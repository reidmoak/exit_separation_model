#!/usr/bin/env python3

# TODO: Figure out equation for tracking from 5000 to 3500 feet (wing with a
#       constant angle with respect to the ground with an initial vert/horiz
#       speed based on whatever it is at z = 5000)
# TODO: Add more configurability with respect to group types (ability to change
#       number of belly fliers, freefliers, etc.)
# TODO: Verify/fix freefly CD and A values to reflect real data
# TODO: Add more granular winds aloft layers to reflect different wind speeds
#       at different altitudes
# TODO: Clean up / Comment code
# TODO: Add user interface to allow for run-time changes to configurable values
# TODO: Add option to import data from CSV file of the load, including number
#       of groups, discipline for each group, average mass of the group, etc.

import time
import sys
import signal
import os
import numpy as np
import matplotlib.pyplot as plt
import subprocess as sp

import params
import const
from skydiver import Skydiver

def signal_handler(sig, frame):
    print('\nYou pressed Ctrl+C, killing all display processes.')
    os.system("killall -s 9 eog")
    sys.exit(0)


# SPECIAL PARAMETERS
rho = params.rho # Air density (NOTE: will eventually not be constant)

# JUMP PARAMETERS
EXIT_ALT     = params.EXIT_ALT      # Exit altitude in feet
BREAKOFF_ALT = params.EXIT_ALT      # Breakoff altitude in feet
PULL_ALT     = params.PULL_ALT      # Pull altitude in feet
IDEAL_SEP    = params.IDEAL_SEP     # Ideal exit separation in feet
SIM_TIME     = params.SIM_TIME      # Number of seconds to run simulation
z0           = params.z0            # Exit altitude (13,500 ft) converted to meters
m            = params.m             # Jumper exit weight (mass, in kg)
Va           = params.Va            # Jump run airspeed in m/s
V_upper      = params.V_upper       # Uppers airspeed in m/s @TODO Make this altitude-dependent
t_sep        = params.t_sep         # Exit separation in seconds
num_groups   = params.num_groups    # Number of groups on the load

# Global runtime variables
Q = 0                               # rho*A*CD, variable to keep code clean
Vg = Va - V_upper                   # Aircraft groundspeed in m/s
x_offset = t_sep * Vg               # Exit separation in meters


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
    plt.axhline(y=BREAKOFF_ALT, color='y', linestyle='dotted', alpha=0.7)
    plt.axhline(y=PULL_ALT, color='r', linestyle='--', alpha=0.7)
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
        pull_x_0_idx = find_nearest(trajectories[i]['z'], PULL_ALT)
        pull_x_1_idx = find_nearest(trajectories[i+1]['z'], PULL_ALT)
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
    if min_sep_dist < IDEAL_SEP:
        plt.hlines(PULL_ALT+200, min_dist_0, min_dist_1, colors='r', label=str(min_sep_dist), linestyle='--')
        plt.text((min_dist_0 + min_dist_1)/2 - const.BLAZE_IT, PULL_ALT+300, str(round(min_sep_dist)) + " ft", color='r')
    else:
        plt.hlines(PULL_ALT+200, min_dist_0, min_dist_1, colors='g', label=str(min_sep_dist), linestyle='--')
        plt.text((min_dist_0 + min_dist_1)/2 - const.BLAZE_IT, PULL_ALT+300, str(round(min_sep_dist)) + " ft", color='g')

    if max_sep_dist < IDEAL_SEP:
        plt.hlines(PULL_ALT+200, max_dist_0, max_dist_1, colors='r', label=str(max_sep_dist), linestyle='--')
        plt.text((max_dist_0 + max_dist_1)/2 - const.BLAZE_IT, PULL_ALT+300, str(round(max_sep_dist)) + " ft", color='r')
    else:
        plt.hlines(PULL_ALT+200, max_dist_0, max_dist_1, colors='g', label=str(max_sep_dist), linestyle='--')
        plt.text((max_dist_0 + max_dist_1)/2 - const.BLAZE_IT, PULL_ALT+300, str(round(max_sep_dist)) + " ft", color='g')

    # Save figure to file and return 0
    plt.savefig("z_vs_x_all_groups.png")
    return 0

if __name__ == "__main__":
    # Initialize signal handler
    signal.signal(signal.SIGINT, signal_handler)

    # Array of time stamps for jump data in seconds
    t = np.array(range(SIM_TIME))

    # Create an array of Skydivers based on num_groups
    load = []
    for i in range(num_groups):
        load.append(Skydiver(m))

    # Compute x and z positions and velocities for each group on the load
    trajectories = []
    for i, jumper in enumerate(load):
        if i < 3: # TODO: Change this to actually function correctly without being hardcoded
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
        u = u - V_upper
        x = x - V_upper*t

        trajectories.append({
            "u" : u*const.M_TO_FT,
            "w" : w*const.M_TO_FT,
            "x" : x*const.M_TO_FT,
            "z" : z*const.M_TO_FT
        })
        
    # print(trajectories)

    DEBUG = 0
    if DEBUG == 1:
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
