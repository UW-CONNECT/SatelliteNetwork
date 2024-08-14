'''
doppler_sim :: Jackie Schellberg :: Last Updated 12/15/2023 

For generating doppler curves of a satellite from provided TLE 

Creates 100 samples/sec, which can be interpolated and used in GNURADIO. 
'''

from skyfield.api import EarthSatellite, Topos,load
import skyfield.api
import numpy as np
import pytz
from datetime import datetime
from dateutil.relativedelta import relativedelta
from matplotlib import pyplot as plt
import pickle

outfile = "ISS_doppler_sim_437_8.pkl"
time_next_pass = 23.9*60
t_end_of_pass = 20
sat = EarthSatellite('1 25544U 98067A   23349.38384836  .00022434  00000+0  39396-3 0  9991',
                     '2 25544  51.6414 151.0242 0002072  34.3439  48.5014 15.50428682429890')
time_scale = skyfield.api.load.timescale()
signal_freq = 437.8e6 

ts = load.timescale()
tz = pytz.timezone('UTC')
dt = tz.localize(datetime.utcnow()+ relativedelta(minutes=time_next_pass))
t0 = ts.utc(dt)
t1 = ts.utc(dt + relativedelta(minutes=t_end_of_pass))      
#npts = t*60*200000  
npts = t_end_of_pass*60*100  
time0 = ts.linspace(t0, t1, npts)

#my_loc = Topos('39.0 N', '105.0 W')
my_loc = Topos('43.0722 N', '89.4008 W')
#my_loc = [43.0722, 89.4008]

def get_radial_velocity(tx_pos, rx_pos):
    ''' see https://github.com/skyfielders/python-skyfield/issues/590 '''
    v_rad = []
    for (tx,rx) in zip(tx_pos,rx_pos): 
        # direct line of sight
        slant_icrf  = tx - rx
        slant       = slant_icrf.position.km
        unit_slant  = slant / np.linalg.norm(slant)

        # satellite velocity
        vel0        = tx.velocity.km_per_s
        v_rel0      = len(vel0)
        unit_vel0   = vel0 / np.linalg.norm(vel0)

        # ground station velocity
        vel1        = rx.velocity.km_per_s
        v_rel1      = len(vel1)
        unit_vel1   = vel1 / np.linalg.norm(vel1)

        # angular deviation of the satellite from the direct line of sight
        theta0  = np.arccos(np.dot(unit_vel0, unit_slant))

        # angular deviation of the ground station from the direct line of sight
        theta1  = np.arccos(np.dot(unit_vel1, unit_slant))

        # radial velocity
        #print(v_rel0 * np.cos(theta0) + v_rel1 * np.cos(theta1))
        v_rad.append( (v_rel0 * np.cos(theta0) + v_rel1 * np.cos(theta1)) * 1000)

    return v_rad

def get_doppler(v_rad, signal_freq):
    doppler = []
    emw = 3e8  
    for v in v_rad: 
        doppler.append( (emw / (emw + v))*signal_freq   )
    return doppler 
    
    
#print(sat.at(time0)[0].velocity)
radial_v = get_radial_velocity(sat.at(time0), my_loc.at(time0))
doppler = np.array(get_doppler(radial_v, signal_freq)) - signal_freq

#print(time0) 

with open(outfile, 'wb') as f:  # Python 3: open(..., 'wb')
    pickle.dump([time0, doppler, t_end_of_pass, signal_freq], f)
   
# doppler values need to be continuous to beu sed 
doppler = np.concatenate((doppler, np.flip(doppler)))
   
print(max(doppler))
plt.figure(1)
xpts = range(0,len(doppler))
plt.plot(xpts, doppler)
plt.show()  


