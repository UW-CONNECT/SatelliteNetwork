from matplotlib import pyplot as plt
import pickle
import numpy as np 
import math
#with open('objs.pkl','rb') as f:  # Python 3: open(..., 'rb')
f = open('ISS_doppler_sim_437_8.pkl','rb')
time0, doppler, t_end_of_pass, signal_freq = pickle.load(f)
f.close()
'''
plt.figure(1)
xpts = range(0,len(doppler))
plt.plot(xpts, doppler)
'''
 
doppler = doppler -signal_freq

fsdop = 100
t = np.arange(0, len(doppler)) / (fsdop)
t = np.arange(0, len(doppler)) / 200000
'''
plt.figure(2)
xpts = range(0,len(doppler))
plt.plot(t, doppler)
plt.show()  
'''


dopplert = 50000
doppler_shift = np.exp(1j * 2 * math.pi * dopplert * t) 
plt.figure(2)
plt.specgram(doppler_shift,Fs=200000)
#plt.plot(np.real(doppler_shift))
#print(len(spec))
plt.show()