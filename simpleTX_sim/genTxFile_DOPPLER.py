'''
tx_stream :: Jackie Schellberg :: Last updated Nov 1, 2023 

For streaming a single set of 'symbols' over UDP to GNUradio 
'''

import socket 
from css_mod import CssMod
from matplotlib import pyplot as plt
import numpy as np
import time
import math
import zmq
import random
import pickle
from scipy.interpolate import CubicSpline
from scipy.interpolate import pchip_interpolate
from scipy import signal
#with open('objs.pkl','rb') as f:  # Python 3: open(..., 'rb')
f = open('doppler_sim_scripts/objs.pkl','rb')
time0, doppler, t_end_of_pass, signal_freq = pickle.load(f)
f.close()


#SF = 7 
SF = 9
N = 2**SF
UPSAMP = 10;
FS = 200000
#symbols = np.concatenate((np.ones(8), range(1,50)))
#symbols =  range(1,50)
#symbols = random.sample(range(1, N), 1000)
'''
plt.figure(1)
xpts = range(0,len(doppler))
plt.plot(xpts, doppler)
plt.show()
'''
doppler_zen = int(np.argwhere(doppler == min(abs(doppler ))) *2000)

fileout_name = "doppler_tests/SF9_100s"

symbols = []
#for ss in range(0,1000):
#    symbols.append(random.randint(1,N-1))
#symbols = np.ones(1000)*5
symbols = np.ones(100)*5

#preamble = [0,0,0,0,0,0,0,0]
#end_delimeter = [0,0,0,0,0,0,0,0] 
#preamble = [1, 0,0,1,1,0,0,1]
preamble = [1,1]
#preamble = [1,1,1,1]
#end_delimeter = [1,0,5,0,0,10,1,1] 
end_delimeter = [3,3,3,3] 
css_modulator = CssMod(N, UPSAMP, preamble, end_delimeter) 
output = css_modulator.symbol2packet(symbols)

# Add doppler 
t = np.linspace(0, len(output)/FS, len(output))
freq_shift = doppler #- signal_freq

#print(len(freq_shift))
#print(len(freq_shift)/100 * FS)
#freq_shift = signal.resample(freq_shift, int(len(freq_shift)/100 * FS))
#freq_shift = signal.resample(freq_shift, int(len(freq_shift)/100 * FS))
print("Expected size: ", int(len(freq_shift)/100 * FS)) 
freq_shift = pchip_interpolate(range(0,(int(len(freq_shift)/100 * FS)), 2000),freq_shift,range(0,(int(len(freq_shift)/100 * FS)))) 
print("ACTUAL: ", len(freq_shift))
#freq_shift = cs(range(0,(int(len(freq_shift)/100 * FS))) ,3)

plt.figure(1)
freq_shift_cut = freq_shift[::(N*10)]
print(len(freq_shift_cut), len(freq_shift))
diffpts = abs(np.diff(freq_shift_cut))
plt.plot(diffpts)

print("MAX DIFF: ", max(diffpts))
print("MIN DIFF: ", min(diffpts))
plt.title("CUT")

plt.figure(2)
plt.plot(freq_shift)
plt.show()  
#print()


#print(len(freq_shift))
#freq_shift = freq_shift[:len(t)]
print(doppler_zen)
halfl = int(len(t)/2)
freq_shift = freq_shift[doppler_zen-halfl:doppler_zen+halfl  ]
print(max(freq_shift))
print("LENGTHS: ", len(freq_shift), len(t))
plt.figure(1)
xpts = range(0,len(freq_shift))
plt.plot(xpts, freq_shift)
plt.show()  
#freq_shift = freq_shift[:len(t)]
output = output * np.exp(1j * 2 * math.pi * freq_shift * t)
plt.figure(3)
plt.specgram(output)
plt.show()

bin_dat = np.float32(css_modulator.ang2bin(output))
#print(bin_dat)


bsl = len(bin_dat)
#bin_dat = np.concatenate((bin_dat,bin_dat,bin_dat,bin_dat,bin_dat))

#print(bin_dat.shape)
#bin_dat = np.concatenate((np.zeros((bsl)), bin_dat))
#print(bin_dat.shape)
#print(len(bin_dat))
bin_dat = np.append(np.float32(np.zeros((bsl))),bin_dat)
bin_dat = bin_dat.tobytes()

#print(np.concatenate((preamble, [49],np.squeeze(symbols), end_delimeter)))

file = open(fileout_name, 'bw')
file.write(bin_dat)
file.close()

file = open(fileout_name + '_ref','w')
file.write(f"{symbols}")
file.close()
