#!/usr/bin/env python3

from urllib.request import urlopen
import numpy as np
import const

def get_ACY_forecast(url):
    for line in file:
        decoded_line = line.decode("utf-8")
        if "ACY " in decoded_line:
            return decoded_line

url = "https://www.aviationweather.gov/windtemp/data?level=low&fcst=06&region=bos&layout=on&date="
file = urlopen(url)

data_points = get_ACY_forecast(file).strip().split(' ')

# Get 3000, 6000, 9000, 12000, and 18000 feet forecasts
altitudes = ['3000', '6000', '9000', '12000', '18000']
forecast = {}
for i, val in enumerate(data_points):
    if i < len(altitudes): forecast[altitudes[i]] = data_points[i+1]
    else: break

winds = {}
for altitude in forecast.keys():
    direction = int(forecast[altitude][0:2])*10
    speed = int(forecast[altitude][2:4])
    winds[altitude] = (direction, speed)

def get_forecast():
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
            exit(69)
    return adj_winds

# Return tuple of "x" and "y" component of winds aloft, relative to jump run
# reference frame -- that is, let jump_run angle be 0 degrees, give relative
# u and v wind velocities
def compute_wind_adj(jump_run, altitude, winds):
    # Adjust winds to be in same reference as jump run angle
    winds = compass_rotate(winds)

    # Change altitude to feet (TODO: how to avoid losing precision?)
    altitude = altitude * const.M_TO_FT    

    # Interpolate linearly between closest angles to get theta for altitude
    if altitude >= 18000:
        theta = np.deg2rad(np.abs(winds['18000'][0] - jump_run))
        speed = winds['18000'][1]
    elif altitude < 18000:
        theta_w = np.interp(altitude, [12000, 18000], \
            [winds['12000'][0], winds['18000'][0]])
        theta = np.deg2rad(np.abs(theta_w - jump_run))
        speed = np.interp(altitude, [12000, 18000], \
            [winds['12000'][1], winds['18000'][1]])
    elif altitude == 12000:
        theta = np.deg2rad(np.abs(winds['12000'][0] - jump_run))
        speed = winds['12000'][1]
    elif altitude < 12000:
        theta_w = np.interp(altitude, [9000, 12000], \
            [winds['9000'][0], winds['12000'][0]])
        theta = np.deg2rad(np.abs(theta_w - jump_run))
        speed = np.interp(altitude, [9000, 12000], \
            [winds['9000'][1], winds['12000'][1]])
    elif altitude == 9000:
        theta = np.deg2rad(np.abs(winds['9000'][0] - jump_run))
        speed = winds['9000'][1]
    elif altitude < 9000:
        theta_w = np.interp(altitude, [6000, 9000], \
            [winds['6000'][0], winds['9000'][0]])
        theta = np.deg2rad(np.abs(theta_w - jump_run))
        speed = np.interp(altitude, [6000, 9000], \
            [winds['6000'][1], winds['9000'][1]])
    elif altitude == 6000: # Don't want to lose precision from interpolating...
        theta = np.deg2rad(np.abs(winds['6000'][0] - jump_run))
        speed = winds['6000'][1]
    elif altitude < 6000:
        theta_w = np.interp(altitude, [3000, 6000], \
            [winds['3000'][0], winds['6000'][0]])
        theta = np.deg2rad(np.abs(theta_w - jump_run))
        speed = np.interp(altitude, [3000, 6000], \
            [winds['3000'][1], winds['6000'][1]])
    elif altitude <= 3000:
        theta = np.deg2rad(np.abs(winds['3000'][0] - jump_run))
        speed = winds['3000'][1]

    #print("altitude = " + str(altitude) + ", theta = " + str(np.rad2deg(theta)) + ", speed = " + str(speed))

    speed = speed * const.KT_TO_MPS
        
    return (speed*np.cos(theta), speed*np.sin(theta))

def print_winds(winds):
    degree_sign = u"\N{DEGREE SIGN}"
    print("\nWinds Aloft for ACY:\n")
    for altitude in winds:
        print(str(altitude) + " ft: " + str(winds[altitude][1])+ " kts from " \
                + str(winds[altitude][0]) + degree_sign)

