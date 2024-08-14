'''
track_multipleTx.py :: Jackie Schellberg :: November 10, 2023 

For tracking connectivity of multiple ground stations 

'''
from SatTracker import SatTracker 
from matplotlib import pyplot as plt
import matplotlib
import numpy as np

# names of the satellites of interest 
# satnames = ['NOAA 15','NOAA 18', 'NOAA 19']

''' Make sure to check the ham radio page if the transponder is active/inactive '''
# satnames = ['NOAA 15','NOAA 18', 'NOAA 19']
satnames = ['FUNCUBE-1 (AO-73)','JY1SAT (JO-97)', 'OSCAR 7 (AO-7)', 'RS-44 & BREEZE-KM R/B', 'XW-2B']
#satnames = ['NOAA 18']

# location of our ground station, longitude and latitude 
# Note that these coordinates suite WGS
madison = [43.0722, 89.4008]
georgia = [33.753746, 84.386330] 

# time in minutes for tracking 
duration = 60*24

granularity = 1
npts = int(duration/granularity) 

# deg above horizon we consider acceptable link location 
minElevationAngle = 25

satTracker = SatTracker(satnames, minElevationAngle)
access, timescales = satTracker.test_access_GS_SAT_GS(madison, georgia, duration, npts)

#ts = np.linspace(0, duration/60, npts)
ts = np.linspace(0, duration/60, npts)

# plot mutual access 
for sat_access in access: 
    plt.figure(1)
    # plt.plot(ts, sat_access)
    plt.scatter(ts[np.argwhere(sat_access > 0)], sat_access[np.argwhere(sat_access > 0)])
plt.legend(satnames)
plt.xlabel("Hours from Current Time")
plt.ylabel("Pass")
plt.title("Access to Madison and Georgia")

access,timescales = satTracker.test_access_GS_SAT(madison, duration, npts)
for sat_access in access: 
    plt.figure(0)
    # plt.plot(ts, sat_access)
    plt.scatter(ts[np.argwhere(sat_access > 0)], sat_access[np.argwhere(sat_access > 0)])
plt.legend(satnames)
plt.xlabel("Hours from Current Time")
plt.ylabel("Pass")
plt.title("Access to Madison only")
plt.show()