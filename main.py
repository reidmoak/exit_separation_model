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

import time
import sys
import signal
import os
import numpy as np
import matplotlib.pyplot as plt
import subprocess as sp
from skydiver import Skydiver

def signal_handler(sig, frame):
    print('\nYou pressed Ctrl+C, killing all display processes.')
    os.system("killall -s 9 eog")
    sys.exit(0)

# CONVERSION CONSTANTS
KT_TO_MPS = 0.51444444444444
LB_TO_KG = 0.45359237
M_TO_FT = 3.280839895
MPS_TO_MPH = 2.2369362920544025

# CONSTANTS - THESE SHOULD LIKELY STAY THE SAME
rho = 1.0               # Air density (will eventually not be constant)
A_HD = 0.35             # Head Down surface area (see http://labman.phys.utk.edu/phys221core/modules/m3/drag.html) 
A_BELLY = 16.0*68.0 / 144.0 / ((M_TO_FT)**2) * 0.7 # Human body surface area -- Estimated by taking my width*height in inches, 
                                                                               #converting to feet, then converting to meters 
                                                                               #and multiplying by an estimate of the percentage 
                                                                               #of the rectangular surface that would be filled 
                                                                               #by a human body falling with an arch
CD_HD = 0.8             # Drag Coefficient of a Head Down flyer (see https://tinyurl.com/ya99r6b4)
CD_BELLY = 1.0          # Estimated drag coefficient of a belly flyer
g = 9.81                # Gravitational acceleration constant in m/s^2
z0 = 13500 * 1/M_TO_FT  # Initial height (13,500 ft) converted to meters

################################################################################
# CONFIGURABLE CONSTANT PARAMETERS #############################################
################################################################################
m = 160*LB_TO_KG        # Jumper exit weight (mass), converted from lbs to kg
Va = 70*KT_TO_MPS       # Aircraft airspeed, converted from knots to m/s
V_upper = 20*KT_TO_MPS  # Average winds aloft speed, converted from knots to m/s
t_sep = 10               # Exit separation in seconds
num_groups = 7          # Number of groups on the load
################################################################################
################################################################################
################################################################################

# Global runtime variables
Q = 0                   # rho*A*CD, variable to keep code clean
Vg = Va - V_upper       # Aircraft groundspeed in m/s
x_offset = t_sep * Vg   # Exit separation in meters

def find_nearest(array, value):
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

    # Full load, all groups
    plt.figure(1)
    for i, jumper in enumerate(trajectories):
        plt.plot(jumper['x'], jumper['z'])
    plt.axhline(y=5000, color='y', linestyle='dotted', alpha=0.7)
    plt.axhline(y=3500, color='r', linestyle='--', alpha=0.7)
    plt.axhline(y=0, color='k')
    plt.grid(alpha=0.4,linestyle='--')
    plt.title("Jump Run Trajectories")
    plt.xlabel("Jump Run (ft)")
    plt.ylabel("Altitude (ft)")

    min_sep_dist = 99999
    min_dist_0 = 0
    min_dist_1 = 0
    max_sep_dist = -99999
    max_dist_0 = 0
    max_dist_1 = 0

    for i, _ in enumerate(trajectories):
        pull_x_0_idx = find_nearest(trajectories[i]['z'], 3500)
        pull_x_1_idx = find_nearest(trajectories[i+1]['z'], 3500)
        pull_x_0 = trajectories[i]['x'][pull_x_0_idx]
        pull_x_1 = trajectories[i+1]['x'][pull_x_1_idx]
        distance = pull_x_1 - pull_x_0
        if round(distance, 4) < round(min_sep_dist, 4):
            min_sep_dist = distance
            min_dist_0 = pull_x_0
            min_dist_1 = pull_x_1
        elif round(distance, 4) > round(max_sep_dist, 4):
            max_sep_dist = distance
            max_dist_0 = pull_x_0
            max_dist_1 = pull_x_1
             
        if i == len(trajectories) - 2:
            break

    if min_sep_dist < 1000:
        plt.hlines(3700, min_dist_0, min_dist_1, colors='r', label=str(min_sep_dist), linestyle='--')
        plt.text((min_dist_0 + min_dist_1)/2 - 420, 3800, str(round(min_sep_dist)) + " ft", color='r')
    else:
        plt.hlines(3700, min_dist_0, min_dist_1, colors='g', label=str(min_sep_dist), linestyle='--')
        plt.text((min_dist_0 + min_dist_1)/2 - 420, 3800, str(round(min_sep_dist)) + " ft", color='g')

    if max_sep_dist < 1000:
        plt.hlines(3700, max_dist_0, max_dist_1, colors='r', label=str(max_sep_dist), linestyle='--')
        plt.text((max_dist_0 + max_dist_1)/2 - 420, 3800, str(round(max_sep_dist)) + " ft", color='r')
    else:
        plt.hlines(3700, max_dist_0, max_dist_1, colors='g', label=str(max_sep_dist), linestyle='--')
        plt.text((max_dist_0 + max_dist_1)/2 - 420, 3800, str(round(max_sep_dist)) + " ft", color='g')

    plt.savefig("z_vs_x_all_groups.png")

    return 0

if __name__ == "__main__":
    # Initialize signal handler
    signal.signal(signal.SIGINT, signal_handler)

    # Array of time stamps for jump data in seconds
    t = np.array(range(80))

    # Create an array of Skydivers based on num_groups
    load = []
    for i in range(num_groups):
        load.append(Skydiver(m))

    # Compute x and z positions and velocities for each group on the load
    trajectories = []
    for i, jumper in enumerate(load):
        if i < 3: # TODO: Change this to actually function correctly without being hardcoded
            A = A_BELLY
            CD = CD_BELLY
            Q = rho*A*CD
            u = jumper.compute_u(t, Va, Q)
            w = jumper.compute_w(t, Va, Q)
            x = jumper.compute_x(t, i*x_offset, Va, Q)
            z = jumper.compute_z(t, z0, Va, Q)
        else:
            A = A_HD
            CD = CD_HD
            Q = rho*A*CD
            u = jumper.compute_u(t, Va, Q)
            w = jumper.compute_w(t, Va, Q)
            x = jumper.compute_x(t, i*x_offset, Va, Q)
            z = jumper.compute_z(t, z0, Va, Q)

        # Account for drift due to uppers    
        u = u - V_upper
        x = x - V_upper*t

        trajectories.append({
            "u" : u*M_TO_FT,
            "w" : w*M_TO_FT,
            "x" : x*M_TO_FT,
            "z" : z*M_TO_FT
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
    plt.plot(t, trajectories[0]['w']*MPS_TO_MPH)
    plt.grid(alpha=.4,linestyle='--')
    plt.title("Vertical Speed vs. Time (First Group)")
    plt.xlabel("Time (s)")
    plt.ylabel("Vertical Speed (mph)")
    plt.savefig("vert_speed_first_grp.png")

    # Single jumper horizontal speed
    plt.figure(3)
    plt.plot(t, trajectories[0]['u']*MPS_TO_MPH)
    plt.grid(alpha=.4,linestyle='--')
    plt.title("Horizontal Speed vs. Time (First Group)")
    plt.xlabel("Time (s)")
    plt.ylabel("Horizontal Speed (mph)")
    plt.savefig("horiz_speed_first_grp.png")

    # Open one of the plots with eog (use arrows to see other PNGs)
    processes.append(sp.Popen("eog z_vs_x_all_groups.png &", shell=True))

    signal.pause()
