#!/usr/bin/env python3

import const

EXIT_ALT        = 13500                     # Exit altitude in feet
BREAKOFF_ALT    = 5000                      # Breakoff altitude in feet
PULL_ALT        = 3500                      # Pull altitude in feet
IDEAL_SEP       = 1000                      # Ideal exit separation (1000 feet)

SIM_TIME        = 120                       # Simulation time in seconds

z0              = EXIT_ALT / const.M_TO_FT  # Initial Altitude in meters

weight          = 160                       # Jumper exit weight in pounds
m               = weight * const.LB_TO_KG   # Jumper mass in kg

rho             = 1.0                       # Air density (eventually not const)

V_upper         = 10 * const.KT_TO_MPS      # Uppers in knots, converted to m/s
Va              = 70 * const.KT_TO_MPS      # Jump run airspeed in knots, 
                                            # converted to m/s

t_sep           = 10                        # Exit separation in seconds
num_groups      = 7                         # Number of groups on the load
