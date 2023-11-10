'''
track_singleTxRx.py :: Jackie Schellberg :: November 9, 2023 

For tracking connectivity between a single ground station and satellites. 

'''
from SatTracker import SatTracker 
from matplotlib import pyplot as plt
import matplotlib
import numpy as np

# names of the satellites of interest 
satnames = ['NOAA 15','NOAA 18', 'NOAA 19']

# location of our ground station, longitude and latitude 
# Note that these coordinates suite WGS
txlatlong = [43.0722, 89.4008]

# time in minutes for tracking 
duration = 60*24

granularity = 1
npts = int(duration/granularity) 

# deg above horizon we consider acceptable link location 
minElevationAngle = 25

satTracker = SatTracker(satnames, minElevationAngle)
access, timescales = satTracker.test_access_GS_SAT(txlatlong, duration, npts)

#ts = np.linspace(0, duration/60, npts)
ts = np.linspace(0, duration/60, npts)

for sat_access in access: 
    plt.figure(1)
    plt.plot(ts, sat_access)
plt.legend(satnames)
plt.xlabel("Hours from Current Time")
plt.ylabel("Pass")
plt.show()
