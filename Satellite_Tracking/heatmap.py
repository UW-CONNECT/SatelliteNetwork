'''
heatmap.py :: Jackie Schellberg :: November 10, 2023 

Replicating the Matlab implementation of heatmap.m : 
Creates a circle about a given lat/long location and shows relative connectivity 
'''
from SatTracker import SatTracker 
from matplotlib import pyplot as plt
import matplotlib
import numpy as np
import pandas as pd 
from mpl_toolkits.basemap import Basemap
import math 
from itertools import chain

# for drawing on the basemap 
def draw_map(m, scale=0.2):
    '''
    https://jakevdp.github.io/PythonDataScienceHandbook/04.13-geographic-data-with-basemap.html
    '''   
    # draw a shaded-relief image
    m.shadedrelief(scale=scale)
    
    # lats and longs are returned as a dictionary
    lats = m.drawparallels(np.linspace(-90, 90, 13))
    lons = m.drawmeridians(np.linspace(-180, 180, 13))

    # keys contain the plt.Line2D instances
    lat_lines = chain(*(tup[1][0] for tup in lats.items()))
    lon_lines = chain(*(tup[1][0] for tup in lons.items()))
    all_lines = chain(lat_lines, lon_lines)
    
    # cycle through these lines and set the desired style
    for line in all_lines:
        line.set(linestyle='-', alpha=0.3, color='w')
    
# names of the satellites of interest (only one for now)
# satnames = ['NOAA 15']
satnames = ['FUNCUBE-1 (AO-73)','JY1SAT (JO-97)', 'OSCAR 7 (AO-7)', 'RS-44 & BREEZE-KM R/B', 'XW-2B']

# location of our ground station, longitude and latitude 
# Note that these coordinates suite WGS
txlatlong = [43.0722, 89.4008]
#txlatlong = [36.6002,121.8947]

# time in minutes for tracking 
duration = 60*24

granularity = 1
npts = int(duration/granularity) 

# deg above horizon we consider acceptable link location 
minElevationAngle = 25

satTracker = SatTracker(satnames, minElevationAngle)
# atlanta, charlotte, nyc 
#rxlatlong = [[33.753746, 84.386330], [35.2271,80.8431] ,[40.7128,74.0060],[51.5072, 0.1276]]

fifty_km_minutes = 0.4522 *2# taken from heatmap.m; not sure what this is 
#fifty_km_minutes = 50
latlongpts = satTracker.getCirclePoints(np.array(txlatlong), 16, fifty_km_minutes)

''' Now test the access, and set coloring according to total minutes available '''
color_list = [] 
# TODO: this only works for a single satellite
for latlong in latlongpts: 
    # test access at a point, and sum access time 
    access, timescales = satTracker.test_access_GS_SAT(latlong, duration, npts)
    #print(len(access))
    access_time = np.sum(access)
    #print(access_time) 
    color_list.append(access_time)

m = Basemap(projection='cyl', resolution=None,
            llcrnrlat=15, urcrnrlat=90,
            llcrnrlon=-140, urcrnrlon=-55, )
draw_map(m)
plt.scatter(latlongpts[:,1], latlongpts[:,0], c=color_list)
plt.colorbar(label="Minutes of Connectivity")
plt.show()