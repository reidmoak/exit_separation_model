#!/usr/bin/env python3

# TODO: Figure out equation for tracking from 5000 to 3500 feet (wing with a
#       constant angle with respect to the ground with an initial vert/horiz
#       speed based on whatever it is at z = 5000)
# TODO: Verify/fix freefly CD and A values to reflect real data
#       NOTE: This will only be possible once I have reliable flysight data for
#       a solo belly jump with known parameters (weight, winds, etc.)
# TODO: Verify that I am doing x calculation correctly, using trapezoid rule of
#       u[t] rather than the compute_x function, which only calculates x in the
#       air reference frame
# TODO: Add option to import data from CSV file of the load, including number
#       of groups, discipline for each group, average mass of the group, etc.
# TODO: Make rho a function of altitude -- would need to make a lookup table in 
#       skydiver.py and make it rho[z] with interpolation. Might be difficult,
#       but I guess I can just compute z[t] first, using z[t-1] as the index 
#       for rho[z] -- so z[t] = f(rho[z[t-1]]). Do this on a new branch probs
# TODO: Clean up / Comment code

# Python Built-In imports
import sys
import signal
import os
import numpy as np
import matplotlib.pyplot as plt
import subprocess as sp
from termcolor import colored

# Custom Module imports
import const
import winds as wind
from parameters import Params
from skydiver import Skydiver


# Special Parameters
rho = 18.436 * (10**-4) # Air density in slugs/ft^3

def signal_handler(sig, frame):
    print('\nYou pressed Ctrl+C, killing all display processes.')
    os.system("killall -s 9 eog")
    sys.exit(0)

# Return array index where array value is closest to provided value
def find_nearest_idx(array, value):
    closest = 0
    found = False
    for i, val in enumerate(array):
        if val > value:
            continue
        elif found is False:
            closest = i
            found = True
            error = value - array[closest]
        if np.abs((array[i-1] - value)) < np.abs(error):
            closest = i-1
    return closest

def adjust_for_uppers(u, x, z, jump_run, winds, t_step):
    # Create empty v[t] and y[t] numpy arrays
    v = np.array(range(len(u)), dtype='f')
    y = np.array(range(len(u)), dtype='f')

    for t in range(len(u)):
        # Get velocity adjustment based on velocity of air column at given t, in ft/s
        u_adj = wind.compute_wind_adj(jump_run, z[t], winds)[0]*const.KT_TO_FPS
        v_adj = wind.compute_wind_adj(jump_run, z[t], winds)[1]*const.KT_TO_FPS

        # Adjust u[t] to be GROUND speed over time, adjusting for movement of 
        # wind column, which is captured by u_adj. No force-driven velocity in
        # the y direction, so v[t] equals just v_adj
        u[t] = u[t] + u_adj
        v[t] = v_adj

        # At time zero, x is 0, so no adjustment needs to be made
        if t == 0:
            continue

        # Calculate GROUND distance over time by using the trapezoidal rule,
        # which is just the area under the ground velocity curve
        # np.trapz below is same as u[t-1]*(1 second) + (u[t] - u[t-1])/2
        # NOTE: Changed to the manual method, since I can use t_step to change
        # the dt increment for the integral
        x[t] = x[t-1] + t_step * (u[t-1] + (u[t] - u[t-1])/2)
        y[t] = y[t-1] + t_step * (v[t-1] + (v[t] - v[t-1])/2)
        #x[t] = x[t-1] + np.trapz([u[t-1], u[t]], [t-1, t], axis=0) 
        #y[t] = y[t-1] + np.trapz([v[t-1], v[t]], [t-1, t], axis=0)

    return u, x, v, y

def plot_y_x(trajectories):
    plt.figure(6)
    for i, jumper in enumerate(trajectories):
        plt.plot(jumper['x'], jumper['y'])
    plt.grid(alpha=0.4,linestyle='--')
    plt.title("Jump Run Trajectories - Bird's Eye View")
    plt.xlabel("Jump Run (ft)")
    plt.ylabel("Perpendicular Drift (ft)")
    plt.savefig("y_vs_x_all_groups.png")

def plot_z_x(trajectories):
    plt.figure(1)
    # Iterate over full load of jumpers, plotting z vs x positions
    for i, jumper in enumerate(trajectories):
        plt.plot(jumper['x'], jumper['z'])
    plt.grid(alpha=0.4,linestyle='--')
    plt.title("Jump Run Trajectories")
    plt.xlabel("Jump Run (ft)")
    plt.ylabel("Altitude (ft)")
    plt.ylim(bottom=0)

    # Plot horizontal lines for breakoff/pull altitudes, as well as the ground
    plt.axhline(y=params.BREAKOFF_ALT, color='y', linestyle='dotted', alpha=0.7)
    plt.axhline(y=params.PULL_ALT, color='r', linestyle='--', alpha=0.7)
    plt.axhline(y=0, color='k')

    # Initialize variables for min and max separation distance
    min_sep_dist = 99999
    min_dist_0 = 0
    min_dist_1 = 0
    max_sep_dist = -99999
    max_dist_0 = 0
    max_dist_1 = 0

    # Iterate over trajectories to find min/max separation distance
    for i, _ in enumerate(trajectories):
        # For each pair of jumpers, search for the x positions corresponding to 
        # the time at which the z position is closest to pull altitude
        pull_x_0_idx = find_nearest_idx(trajectories[i]['z'], params.PULL_ALT)
        pull_x_1_idx = find_nearest_idx(trajectories[i+1]['z'], params.PULL_ALT)
        pull_x_0 = trajectories[i]['x'][pull_x_0_idx]
        pull_x_1 = trajectories[i+1]['x'][pull_x_1_idx]

        # Compute exit separation and update min/max values based on saved ones
        distance = pull_x_1 - pull_x_0
        if round(distance, 4) < round(min_sep_dist, 4):
            min_sep_dist = distance
            min_dist_0 = pull_x_0
            min_dist_1 = pull_x_1
        elif round(distance, 4) > round(max_sep_dist, 4):
            max_sep_dist = distance
            max_dist_0 = pull_x_0
            max_dist_1 = pull_x_1
             
        # Break at num_groups minus 2, to account for 0-based indexing and the
        # fact that we're looking at index i+1
        if i == len(trajectories) - 2:
            break

    # Plot and label horizontal lines for min/max separation distance, using
    # the color red if distance is less than the ideal separation distance
    # (nominally 1000 feet) and green otherwise
    min_color = 'r' if min_sep_dist < params.IDEAL_SEP else 'g'
    plt.hlines(params.PULL_ALT+200, min_dist_0, min_dist_1, colors=min_color, \
               label=str(min_sep_dist), linestyle='--')
    plt.text((min_dist_0 + min_dist_1)/2 - const.BLAZE_IT, params.PULL_ALT+300, \
             str(round(min_sep_dist)) + " ft", color=min_color)
    
    max_color = 'r' if max_sep_dist < params.IDEAL_SEP else 'g'
    plt.hlines(params.PULL_ALT+200, max_dist_0, max_dist_1, colors=max_color, \
               label=str(max_sep_dist), linestyle='--')
    plt.text((max_dist_0 + max_dist_1)/2 - const.BLAZE_IT, params.PULL_ALT+300, \
             str(round(max_sep_dist)) + " ft", color=max_color)

    # Save figure to file and return 0
    plt.savefig("z_vs_x_all_groups.png")
    return 0

# TODO: Figure this one out...
def plot_3d(trajectories):
    plt.figure(7)
    #plt.grid(alpha=0.4,linestyle='--')
    plt.title("Trajectories - 3D")
    x = trajectories[0]['x']
    y = trajectories[0]['y']
    z = trajectories[0]['z']
    ax = plt.axes(projection='3d')
    #plt.xlabel("x (ft)")
    #plt.ylabel("y (ft)")
    #plt.zlabel("z (ft)")
    ax.plot3D(x, y, z, 'gray')
    plt.savefig("3d_plot.png")

def main_menu(winds):
    options = []
    options.append("Print winds aloft (ACY)")
    options.append("Print simulation parameters")
    options.append("Modify simulation parameters")
    options.append("Enable debug options")
    options.append("Set jump run and exit separation to optimal values")
    options.append(colored("Run simulation", 'yellow'))

    rec_jump_rum = None
    rec_t_sep = None

    def print_title():
        os.system('clear')
        print(colored("EXIT SEPARATION MODEL\n", 'green'))

    print_title()

    while(True):
        for i, opt in enumerate(options):
            print(str(i+1) + ") " + opt)
        ans = input("\nEnter selection: ")
        try:
            ans = int(ans)
            if ans < 0 or ans > len(options):
                print_title()
                print(colored("Invalid entry, input must be between 1 and " + \
                              str(len(options)) + ".\n", 'red'))
            elif ans == 1: # Print winds aloft (ACY) 
                print_title()
                rec_jump_run, rec_t_sep = wind.print_winds(winds, params.aircraft, params.EXIT_ALT)
                print("")
            elif ans == 2: # Print simulation paramters
                print_title()
                params.show(False)
            elif ans == 3: # Modify simulation parameters
                params.setup()
                print_title()
            elif ans == 4: # Enable debug options
                print_title()
                print(colored("Option 4 is currently under construction...\n", 'cyan'))
            elif ans == 5: # Set jump run and exit separation to optimal values
                rec_jump_run, rec_t_sep = wind.print_winds(winds, params.aircraft, params.EXIT_ALT)
                params.jump_run = int(rec_jump_run)
                params.t_sep = rec_t_sep
                print_title()
                params.show(False)
            elif ans == 6: # Run simulation
                return
        except ValueError:
            print_title()
            print(colored("Invalid entry, input must be a number.\n", 'red'))
            continue
            
def generic_plot(x, y, fid, title, xlabel, ylabel, filename):
    plt.figure(fid)
    plt.plot(x, y)
    plt.grid(alpha=.4,linestyle='--')
    plt.title(title)
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.savefig(filename)

if __name__ == "__main__":
    # Initialize params class, which will prompt user to setup parameter values
    params = Params()

    # Initialize signal handler
    signal.signal(signal.SIGINT, signal_handler)

    # Winds
    winds = wind.get_forecast()

    # Aircraft speeds TODO: Verify these! Units in knots
    aircraft_speeds = const.AIRCRAFT_SPEEDS

    main_menu(winds)

    # Initialize runtime variables
    Q = 0                                       # rho*A*CD, to keep code clean
    z0 = params.EXIT_ALT                        # Initial Altitude in feet
    m = params.weight / const.g                 # Jumper mass in slugs
    Va = aircraft_speeds.get(params.aircraft)   # Aircraft airspeed in knots
    num_groups = params.num_rw_groups + \
                 params.num_ff_groups
    sim_time = 120                              # TODO: Figure out how to make dynamic 
    t_step = 0.01                               # Seconds between sim time steps

    # Array of time stamps for jump data in seconds
    t = np.array(range(sim_time))
    t = np.array(np.linspace(0, sim_time, int(sim_time/t_step) + 1))

    # Winds Aloft DEBUG - makes them constant at 180 from jump run
    simple_winds = False
    if simple_winds is True:
        V_upper = params.V_upper            # Constant uppers in knots
        Vg = Va - V_upper                   # Aircraft groundspeed in knots
    else:
        V_upper_adj = wind.compute_wind_adj(params.jump_run, params.EXIT_ALT, winds)
        Vg = Va + V_upper_adj[0]            # Aircraft ground speed in knots

    # Compute exit separation distance in feet
    x_offset = params.t_sep * (Vg*const.KT_TO_FPS) 

    # Create an array of Skydivers based on num_groups
    load = []
    for i in range(num_groups):
        load.append(Skydiver(m))

    # Compute x and z positions and velocities for each group on the load,
    # converting speeds from knots to m/s first
    # NOTE: compute_x not really necessary, since we compute the x distance in
    # the ground refernce frame when we call adjust_for_uppers
    trajectories = []
    for i, jumper in enumerate(load):
        if i < params.num_rw_groups:
            A = const.A_BELLY
            CD = const.CD_BELLY
            Q = rho*A*CD
            u = jumper.compute_u(t, Va*const.KT_TO_FPS, Q)
            w = jumper.compute_w(t, Va*const.KT_TO_FPS, Q)
            x = jumper.compute_x(t, i*x_offset, Va*const.KT_TO_FPS, Q)
            z = jumper.compute_z(t, z0, Va*const.KT_TO_FPS, Q)
        else:
            A = const.A_HD
            CD = const.CD_HD
            Q = rho*A*CD
            u = jumper.compute_u(t, Va*const.KT_TO_FPS, Q)
            w = jumper.compute_w(t, Va*const.KT_TO_FPS, Q)
            x = jumper.compute_x(t, i*x_offset, Va*const.KT_TO_FPS, Q)
            z = jumper.compute_z(t, z0, Va*const.KT_TO_FPS, Q)

        # Account for drift due to uppers    
        if simple_winds is True:
            u = u - V_upper*const.KT_TO_FPS
            x = x - V_upper*const.KT_TO_FPS*t # Okay since V_upper is constant
        else:
            u, x, v, y = adjust_for_uppers(u, x, z, params.jump_run, winds, t_step)

        trajectories.append({
            "u" : u*const.FPS_TO_MPH, # Convert to MPH for plots
            "v" : v*const.FPS_TO_MPH, # Convert to MPH for plots
            "w" : w*const.FPS_TO_MPH, # Convert to MPH for plots
            "x" : x,
            "y" : y,
            "z" : z
        })
  
    # Trajectory DEBUG
    traj_debug = False
    if traj_debug is True:
        #print(trajectories)
        for i, z in enumerate(trajectories[0]['z']):
            print("time = " + str(i) + " sec, Altitude = " + str(z) + " ft")
            if z < 0: 
                break
        for i, z in enumerate(trajectories[len(trajectories)-1]['z']):
            print("time = " + str(i) + " sec, Altitude = " + str(z) + " ft")
            if z < 0:
                break
        for i, u in enumerate(trajectories[0]['u']):
            print("time = " + str(i) + " sec, u[t] = " + str(u) + " mph, altitude = " + \
                  str(trajectories[0]['z'][i]) + " ft")

    ###################
    #      PLOTS      #
    ###################
    processes = []

    # Full load, all groups 
    #   z vs x is side view perpendicular to left-to-right jump run 
    #   y vs x is birds-eye view (jump-run still left to right)
    plot_z_x(trajectories)
    plot_y_x(trajectories)

    # Single jumper vertical speed - Belly group
    generic_plot(t, trajectories[0]['w'], 2, \
                 "Vertical Speed vs. Time (First Group)", "Time (s)", \
                 "Vertical Speed (mph)", "vert_speed_first_grp.png")

    # Single jumper horizontal speed - Belly group
    generic_plot(t, trajectories[0]['u'], 3, \
                 "Horizontal Speed vs. Time (First Group)", "Time (s)", \
                 "Horizontal Speed (mph)", "horiz_speed_first_grp.png")

    # Single jumper vertical speed - Freefly group
    generic_plot(t, trajectories[len(trajectories)-1]['w'], 4, \
                 "Vertical Speed vs. Time (Last Group)", "Time (s)", \
                 "Vertical Speed (mph)", "vert_speed_last_grp.png")

    # Single jumper horizontal speed - Freefly group
    generic_plot(t, trajectories[len(trajectories)-1]['u'], 5, \
                 "Horizontal Speed vs. Time (Last Group)", "Time (s)", \
                 "Horizontal Speed (mph)", "horiz_speed_last_grp.png")

    # Open one of the plots with eog (use arrows to see other PNGs)
    processes.append(sp.Popen("eog z_vs_x_all_groups.png &", shell=True))

    signal.pause()

