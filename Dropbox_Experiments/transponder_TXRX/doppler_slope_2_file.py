import numpy as np 
import math 

FS = 200e3 
doppler_slope = 50 # hz/s 
duration = 15*60 # seconds  

t = np.arange(0, duration) / FS 
doppler = np.exp(1j * 2 * math.pi * doppler_slope * t**2/2)




