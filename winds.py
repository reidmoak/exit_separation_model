#!/usr/bin/env python3

from urllib.request import urlopen

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

print(winds)
