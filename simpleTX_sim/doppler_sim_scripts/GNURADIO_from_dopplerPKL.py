import pickle 
from scipy.interpolate import pchip_interpolate
import numpy as np 
import math
from matplotlib import pyplot as plt
import sys
sys.path.insert(0, '..')
from css_mod import CssMod

FS = 200000
FS_dop = 100
infile = "ISS_doppler_sim_437_8.pkl"
nums2 = FS # number of doppler samples (this corresponds to 2 seconds)

# info not needed
SF = 9 
preamble = [] 
end_delimeter = []
N = 2**SF
UPSAMP = 10

f = open(infile,'rb')
time0, doppler, t_end_of_pass, signal_freq = pickle.load(f)
f.close()



print(len(doppler))

doppler = pchip_interpolate(range(0,(int(len(doppler)/FS_dop * FS)), 2000),doppler,range(0,(int(len(doppler)/FS_dop * FS)))) 
print(len(doppler))

midpoint = np.argwhere(np.diff(doppler)==min(np.diff(doppler)))[0]
doppler = doppler[int(midpoint-nums2):int(midpoint+nums2)]
doppler = np.concatenate((doppler,doppler[::-1]))
plt.figure(1)
xpts = range(0,len(doppler))
plt.plot(xpts, doppler)
plt.show()

t = np.linspace(0, len(doppler)/FS, len(doppler))

output = np.exp(1j * 2 * math.pi * doppler * t)
#wsoutput = output[

css_modulator = CssMod(N, UPSAMP, preamble, end_delimeter) 
bin_dat = np.float32(css_modulator.ang2bin_nopad(output))
bin_dat = bin_dat.tobytes()

file = open("GNURADIO_DOP_REF_small", 'bw')
file.write(bin_dat)
file.close()