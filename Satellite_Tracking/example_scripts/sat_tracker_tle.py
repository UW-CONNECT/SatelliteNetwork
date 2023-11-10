'''
This code was from 
https://space.stackexchange.com/questions/58268/satellite-ground-track-calculation 
and updated TLE from https://celestrak.org/NORAD/elements/stations.txt

This is tracking the ISS across a 200 minute timescale

'''

import pytz
from datetime import datetime
from dateutil.relativedelta import relativedelta
from skyfield.api import load, wgs84, EarthSatellite

ts = load.timescale()

# Updated version of the listed two-line elements

line1 = '1 25544U 98067A   23313.48047773  .00055397  00000+0  96054-3 0  9997'
line2 = '2 25544  51.6436 328.8525 0001280  74.6781  77.0191 15.50331811424322'
satellite = EarthSatellite(line1, line2, 'ISS (ZARYA)', ts)
print(satellite)

# Get the current time in a timezone-aware fashion.

tz = pytz.timezone('UTC')
dt = tz.localize(datetime.utcnow())
print(f"Exectution time: {dt:%Y-%m-%d %H:%M:%S %Z}\n")

# Split the next 200 minutes into a collection of 101-evenly
# spaced Timescales (every two minutes, plus endpoints)

t0 = ts.utc(dt)
t1 = ts.utc(dt + relativedelta(minutes=200))
timescales = ts.linspace(t0, t1, 101)

# calculate the subpoints.

geocentrics = satellite.at(timescales)
subpoints = wgs84.subpoint_of(geocentrics)

# Print a nicely-formatted, tab-delimited time, latitude, and longitude.

for t, lat, lon in zip(timescales,
                       subpoints.latitude.degrees,
                       subpoints.longitude.degrees):
    print(f"{t.astimezone(tz):%Y-%m-%d %H:%M:%S}\t{lat:8.2f}\t{lon:8.2f}") 