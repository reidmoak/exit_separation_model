#!/usr/bin/env python3

from urllib.request import urlopen
import numpy as np
import const
from termcolor import colored

# Get winds aloft forecast for ACY
def get_ACY_forecast(url):
    file_lines = urlopen(url)
    for line in file_lines:
        decoded_line = line.decode("utf-8")
        if "ACY " in decoded_line:
            return decoded_line


def get_forecast():
    # Winds aloft URL
    url = "https://www.aviationweather.gov/windtemp/data?level=low&fcst=06&region=bos&layout=on&date="
    altitudes = [3000, 6000, 9000, 12000, 18000]
    winds = {}
    forecast = {}

    data_points = get_ACY_forecast(url).strip().split(' ')

    # Get 3000, 6000, 9000, 12000, and 18000 feet forecasts
    for i, val in enumerate(data_points):
        if i < len(altitudes): forecast[altitudes[i]] = data_points[i+1]
        else: break

    for altitude in forecast.keys():
        direction = int(forecast[altitude][0:2])*10
        speed = int(forecast[altitude][2:4])
        winds[altitude] = (direction, speed)

    return winds

def compass_rotate(winds):
    adj_winds = {}
    for altitude in winds:
        if winds[altitude][0] >= 0 and winds[altitude][0] <=180:
            adj_winds[altitude] = (winds[altitude][0]+180, winds[altitude][1])
        elif winds[altitude][0] > 180 and winds[altitude][0] <=360:
            adj_winds[altitude] = (winds[altitude][0]-180, winds[altitude][1])
        else:
            print("HOW DID I GET HERE???")
            exit(69) # Nice
    return adj_winds 

# Return tuple of "x" and "y" component of winds aloft, relative to jump run
# reference frame -- that is, let jump_run angle be 0 degrees, give relative
# u and v wind velocities
def compute_wind_adj(jump_run, altitude, winds):
    # Adjust winds to be in same reference as jump run angle
    winds = compass_rotate(winds)

    # Interpolate linearly between closest angles to get theta for altitude
    # TODO: Make this less iterative...probably can have a function that looks
    #       for the 2 closest altitudes and then interpolate between?
    if altitude >= 18000:
        theta = np.deg2rad(np.abs(winds[18000][0] - jump_run))
        speed = winds[18000][1]
    elif altitude > 12000:
        theta_w = np.interp(altitude, [12000, 18000], \
            [winds[12000][0], winds[18000][0]])
        theta = np.deg2rad(np.abs(theta_w - jump_run))
        speed = np.interp(altitude, [12000, 18000], \
            [winds[12000][1], winds[18000][1]])
    elif altitude == 12000:
        theta = np.deg2rad(np.abs(winds[12000][0] - jump_run))
        speed = winds[12000][1]
    elif altitude > 9000:
        theta_w = np.interp(altitude, [9000, 12000], \
            [winds[9000][0], winds[12000][0]])
        theta = np.deg2rad(np.abs(theta_w - jump_run))
        speed = np.interp(altitude, [9000, 12000], \
            [winds[9000][1], winds[12000][1]])
    elif altitude == 9000:
        theta = np.deg2rad(np.abs(winds[9000][0] - jump_run))
        speed = winds[9000][1]
    elif altitude > 6000:
        theta_w = np.interp(altitude, [6000, 9000], \
            [winds[6000][0], winds[9000][0]])
        theta = np.deg2rad(np.abs(theta_w - jump_run))
        speed = np.interp(altitude, [6000, 9000], \
            [winds[6000][1], winds[9000][1]])
    elif altitude == 6000: # Don't want to lose precision from interpolating...
        theta = np.deg2rad(np.abs(winds[6000][0] - jump_run))
        speed = winds[6000][1]
    elif altitude > 3000:
        theta_w = np.interp(altitude, [3000, 6000], \
            [winds[3000][0], winds[6000][0]])
        theta = np.deg2rad(np.abs(theta_w - jump_run))
        speed = np.interp(altitude, [3000, 6000], \
            [winds[3000][1], winds[6000][1]])
    elif altitude <= 3000:
        theta = np.deg2rad(np.abs(winds[3000][0] - jump_run))
        speed = winds[3000][1]
       
    # Return u_adj, v_adj in knots
    #print("theta = " + str(np.rad2deg(theta)) + ", altitude = " + str(altitude) + \
    #      ", u_adj = " + str(speed*np.cos(theta)) + ", v_adj = " + \
    #                      str(speed*np.sin(theta)))
    
    return (speed*np.cos(theta), speed*np.sin(theta))

def print_winds(winds, aircraft, exit_alt):
    degree_sign = u"\N{DEGREE SIGN}"

    print("Winds Aloft for ACY:\n")
    for altitude in list(reversed((sorted(winds.keys())))):
        print("\t" + str(altitude) + " ft: " + str(winds[altitude][1])+ " kts from " \
                + str(winds[altitude][0]) + degree_sign)

    # TODO: Make this not hardcoded -- if EXIT_ALT is changed from 13500 to something else,
    #       the code below will not work
    exit_uppers = np.interp(exit_alt, [12000, 18000], [winds[12000][1], winds[18000][1]])
    exit_deg = np.interp(exit_alt, [12000, 18000], [winds[12000][0], winds[18000][0]])

    print("")
    exit_sep = exit_sep_chart(exit_uppers, aircraft)

    print("\tSuggested jump run: " + colored(str(int(exit_deg)) + degree_sign, 'magenta'))

    return exit_deg, exit_sep 

def exit_sep_chart(exit_uppers, aircraft):
    # Using values from https://www.dropzone.com/articles/safety/exit-separation-r876/

    # Get aircraft-dependent TAS (true airspeed), in knots 
    V_tas = const.AIRCRAFT_SPEEDS.get(aircraft)

    # Compute ground speed based on (1) aircraft airspeed and (2) uppers at EXIT_ALT
    V_g = V_tas - exit_uppers
    exit_sep = {
        100: 6,
        95: 7,
        90: 7,
        85: 7,
        80: 8, 
        75: 8,
        70: 9,
        65: 10,
        60: 10,
        55: 11,
        50: 12,
        45: 14,
        40: 15,
        35: 17,
        30: 20,
        25: 24,
        20: 30,
        15: 40,
        10: 60,
        5: 119
    }.get(5 * round(V_g/5)) # Need to round ground speed to nearest mult of 5

    print("\tAircraft ground speed: " + colored(str(V_g) + " kts", 'magenta'))
    print("\tRecommended exit sep: " + colored(str(exit_sep) + " seconds", 'magenta'))

    return exit_sep

