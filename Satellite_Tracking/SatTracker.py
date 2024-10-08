'''
SatTracker.py :: Jackie Schellberg :: 11/9/2023 

References 
https://rhodesmill.org/skyfield/earth-satellites.html 

'''
import pytz
from datetime import datetime
from dateutil.relativedelta import relativedelta
from skyfield.api import load, wgs84, EarthSatellite, N, S, E, W
import numpy as np
import math

class SatTracker: 
    def __init__(self, sats, minElevationAngle):    
        '''
        sats - string names in a list of satellites of interest 
        minElevationAngle - angle in degrees above which we consider connected
        '''
        # Degrees from the horizon that the satellite needs to be from GS for access 
        self.minElevationAngle = minElevationAngle
                
        # load satellite data TLE data for all sats 
        #stations_url = 'http://celestrak.org/NORAD/elements/stations.txt' # this file contains the ISS 
        # stations_url = 'https://celestrak.org/NORAD/elements/gp.php?GROUP=active&FORMAT=tle'    # this file contains NOAA sats
        
        # SATNOGS has most of our satellites of interet. 
        stations_url = 'https://celestrak.org/NORAD/elements/gp.php?GROUP=satnogs&FORMAT=tle'
        # [TODO] If we are interested in SATs across these list, just concat the output 
        # satellites = load.tle_file(stations_url,reload=True)
        # [TODO] Check cached tle_file only if two hours has no yet passed 
        satellites = load.tle_file(stations_url)
        by_name = {sat.name: sat for sat in satellites}
        #satellite = by_name[sats]
        
        self.satellites = [] 
        for sat in sats:
            self.satellites.append(by_name[sat])
                
    def test_access_GS_SAT(self, gs, t, npts):
        '''
        A duration 't' in minutes, return true at a timepoint if gs has LOS to sat according to minElevationAngle
        '''
        ts = load.timescale()
        access_out = [] 
        
        for sat in self.satellites: 
            access_out_curr = [] 
            location = wgs84.latlon(gs[0]* N, gs[1] * W)
            
            tz = pytz.timezone('UTC')
            dt = tz.localize(datetime.utcnow())
            t0 = ts.utc(dt)
            t1 = ts.utc(dt + relativedelta(minutes=t))         
            timescales = ts.linspace(t0, t1, npts)
            
            difference = sat - location 
            topocentric = difference.at(timescales)
            
            # altitude angle and azimuth angle are the same 
            alt, az, distance = topocentric.altaz()
            
            for ti in range(0,len(timescales)): 
                if alt.degrees[ti] > self.minElevationAngle: 
                    access_out_curr.append(1)
                    # access_out_curr.append(alt.degrees[ti])
                else: 
                    access_out_curr.append(0)
                    
            access_out.append(np.squeeze(access_out_curr))
            
        return access_out, timescales 

    def test_access_GS_SAT_GS(self, gs1, gs2, t, npts): 
        '''
        Finds when gs1 and gs2 share simultaneous connectivity with a satellite. 
        '''
        access_out = [] 
        access_gs1, timescales = self.test_access_GS_SAT(gs1, t, npts)
        access_gs2, timescales = self.test_access_GS_SAT(gs2, t, npts) 
        
        for i in range(0, len(access_gs1)):
            tmp_list = np.logical_and(access_gs1[i],access_gs2[i])
            access_out.append(tmp_list)
        
        return access_out, timescales
        #return np.logical_and(access_gs1, access_gs2) 
        
    def test_access_GS_SAT_multiGS(self, gs1, gsn, t, npts):
        '''
        Finds when gs1 (Tx) has access to at least 1 gsn (Rx) 
        '''
        # bad way to init
        access_out,timescales = self.test_access_GS_SAT_GS(gs1, gsn[0], t, npts)
        
        # we have access if ANY rx is connected, so 'or' the results 
        for gs2 in gsn: 
            access,timescales = self.test_access_GS_SAT_GS(gs1, gs2, t, npts)
                        
            for i in range(0, len(access)):
                tmp_list = np.logical_or(access[i],access_out[i])
                access_out.append(tmp_list)
                
        return access_out, timescales
                
    def getCirclePoints(self,center_point, nCirc, delta):
        '''
        Generate points from N, S, E, W, at max distance nCirc * delta (in lat/long) 
        '''
        center_point = np.array([center_point[0]*N, center_point[1] *W])
        points_out = []
        delta_ang = delta/math.sqrt(2)
                
        for i in range(1, nCirc+1):
            # programming horror 
            points_out.append(center_point + [0, i*delta]) 
            points_out.append(center_point + [0, -i*delta]) 
            points_out.append(center_point + [i*delta, 0]) 
            points_out.append(center_point + [-i*delta, 0]) 
            points_out.append(center_point + [delta_ang*i, delta_ang*i])
            points_out.append(center_point + [delta_ang*i, -delta_ang*i])
            points_out.append(center_point + [-delta_ang*i, delta_ang*i])
            points_out.append(center_point + [-delta_ang*i, -delta_ang*i])
        return np.squeeze(points_out) 