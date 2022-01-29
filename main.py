#!/usr/bin/env python3

# TODO: Figure out equation for tracking from 5000 to 3500 feet (wing with a
#       constant angle with respect to the ground with an initial vert/horiz
#       speed based on whatever it is at z = 5000)
# TODO: Verify/fix freefly CD and A values to reflect real data
# TODO: Clean up / Comment code
# TODO: Add option to import data from CSV file of the load, including number
#       of groups, discipline for each group, average mass of the group, etc.
# TODO: Add in y[t] to be able to create 3D plots
# TODO: Add higher up main menu, where you can do other things besides change
#       variables, like print out winds aloft average speed/direction
# TODO: Look into whether I should make rho non-constant -- how much does 
#       air density change with altitude?

# Python Built-In imports
import time
import sys
import signal
import os
import numpy as np
import matplotlib.pyplot as plt
import subprocess as sp
from termcolor import colored
from scipy import integrate

# Custom Module imports
import const
import winds as wind
from parameters import Params
from skydiver import Skydiver


# Special Parameters
rho = 1.0 # Air density (NOTE: will eventually not be constant)

def signal_handler(sig, frame):
    print('\nYou pressed Ctrl+C, killing all display processes.')
    os.system("killall -s 9 eog")
    sys.exit(0)

# Return array index where array value is closest to provided value
def find_nearest(array, value):
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

def adjust_for_uppers(u, x, z, jump_run, winds):
    for t in range(len(u)):
        u_adj = wind.compute_wind_adj(jump_run, z[t], winds)[0]

        # u (x velocity) adjustment DEBUG
        uadj_debug = False
        if uadj_debug is True:
            print("altitude = " + str(z[t]*const.M_TO_FT) + " u[t] = " + str(u[t]) + ", u_adj = " + str(u_adj))
        
        u[t] = u[t] + u_adj
        x[t] = x[t] + u_adj*t 
    return u, x

def plot_trajectories(trajectories):
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
        pull_x_0_idx = find_nearest(trajectories[i]['z'], params.PULL_ALT)
        pull_x_1_idx = find_nearest(trajectories[i+1]['z'], params.PULL_ALT)
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
            elif ans == 1: 
                print_title()
                rec_jump_run, rec_t_sep = wind.print_winds(winds, params.aircraft, params.EXIT_ALT)
                print("")
            elif ans == 2:
                print_title()
                params.show(False)
            elif ans == 3:
                params.setup()
                print_title()
            elif ans == 4:
                print_title()
                print(colored("Option 4 is currently under construction...\n", 'cyan'))
            elif ans == 5:
                if rec_jump_run is None or rec_t_sep is None:
                    rec_jump_run, rec_t_sep = wind.print_winds(winds, params.aircraft, params.EXIT_ALT)
                params.jump_run = rec_jump_run
                params.t_sep = rec_t_sep
                print_title()
                params.show(False)
            elif ans == 6:
                return
        except ValueError:
            print_title()
            print(colored("Invalid entry, input must be a number.\n", 'red'))
            continue
            


if __name__ == "__main__":
    # Initialize params class, which will prompt user to setup parameter values
    params = Params()

    # Initialize signal handler
    signal.signal(signal.SIGINT, signal_handler)

    # Winds
    winds = wind.get_forecast()
    wind.print_winds(winds, params.aircraft, params.EXIT_ALT)

    # Aircraft speeds TODO: Verify these!
    aircraft_speeds = {
        "Caravan"   : 70,
        "Otter"     : 60,
        "Skyvan"    : 90,
        "Cessna 182": 55
    }

    main_menu(winds)

    # Initialize runtime variables
    Q = 0                                       # rho*A*CD, to keep code clean
    z0 = params.EXIT_ALT / const.M_TO_FT        # Initial Altitude in meters
    m = params.weight * const.LB_TO_KG          # Jumper mass in kg
    Va = aircraft_speeds.get(params.aircraft) \
        * const.KT_TO_MPS                       # Aircraft airspeed in m/s   
    num_groups = params.num_rw_groups + params.num_ff_groups
    sim_time = num_groups * 20                  # Simulation time in seconds
    sim_time = 120                              # TODO: Figure out how to make 
                                                # this dynamic, since the above 
                                                # line makes the x limit of the 
                                                # plot very negative

    # Array of time stamps for jump data in seconds
    t = np.array(range(sim_time))

    # Winds Aloft DEBUG - makes them constant at 180 from jump run
    simple_winds = False
    if simple_winds is True:
        V_upper = params.V_upper * const.KT_TO_MPS  # Constat uppers in m/s
        Vg = Va - V_upper # Aircraft groundspeed in m/s
    else:
        V_upper_adj = wind.compute_wind_adj(params.jump_run, \
                                            params.EXIT_ALT/const.M_TO_FT, \
                                            winds)
        Vg = Va + V_upper_adj[0]        
        x_offset = params.t_sep * Vg    # Exit separation in meters

    # Create an array of Skydivers based on num_groups
    load = []
    for i in range(num_groups):
        load.append(Skydiver(m))

    # Compute x and z positions and velocities for each group on the load
    trajectories = []
    for i, jumper in enumerate(load):
        if i < params.num_rw_groups:
            A = const.A_BELLY
            CD = const.CD_BELLY
            Q = rho*A*CD
            u = jumper.compute_u(t, Va, Q)
            w = jumper.compute_w(t, Va, Q)
            x = jumper.compute_x(t, i*x_offset, Va, Q)
            z = jumper.compute_z(t, z0, Va, Q)
        else:
            A = const.A_HD
            CD = const.CD_HD
            Q = rho*A*CD
            u = jumper.compute_u(t, Va, Q)
            w = jumper.compute_w(t, Va, Q)
            x = jumper.compute_x(t, i*x_offset, Va, Q)
            z = jumper.compute_z(t, z0, Va, Q)

        # Account for drift due to uppers    
        if simple_winds is True:
            u = u - V_upper
            x = x - V_upper*t
        else:
            u, x = adjust_for_uppers(u, x, z, params.jump_run, winds)

        trajectories.append({
            "u" : u*const.M_TO_FT,
            "w" : w*const.M_TO_FT,
            "x" : x*const.M_TO_FT,
            "z" : z*const.M_TO_FT
        })
  
    # Trajectory DEBUG
    traj_debug = False
    if traj_debug is True:
        print(trajectories)
        for i, z in enumerate(trajectories[0]['z']):
            print("time = " + str(i) + " sec, Altitude = " + str(z) + " ft")

    ###################
    #      PLOTS      #
    ###################
    processes = []

    # Full load, all groups 
    plot_trajectories(trajectories)

    # Single jumper vertical speed - TODO: Make this its own function (with horiz speed)
    plt.figure(2)
    plt.plot(t, trajectories[0]['w']*const.MPS_TO_MPH)
    plt.grid(alpha=.4,linestyle='--')
    plt.title("Vertical Speed vs. Time (First Group)")
    plt.xlabel("Time (s)")
    plt.ylabel("Vertical Speed (mph)")
    plt.savefig("vert_speed_first_grp.png")

    # Single jumper horizontal speed
    plt.figure(3)
    plt.plot(t, trajectories[0]['u']*const.MPS_TO_MPH)
    plt.grid(alpha=.4,linestyle='--')
    plt.title("Horizontal Speed vs. Time (First Group)")
    plt.xlabel("Time (s)")
    plt.ylabel("Horizontal Speed (mph)")
    plt.savefig("horiz_speed_first_grp.png")

    # Open one of the plots with eog (use arrows to see other PNGs)
    processes.append(sp.Popen("eog z_vs_x_all_groups.png &", shell=True))

    signal.pause()

