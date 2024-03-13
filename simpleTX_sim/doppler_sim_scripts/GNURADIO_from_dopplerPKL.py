import pickle 
from scipy.interpolate import pchip_interpolate
import numpy as np 
import math
from matplotlib import pyplot as plt
import sys
sys.path.insert(0, '..')
from css_mod import CssMod

#FS = 200000
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
BW = 20000

f = open(infile,'rb')
time0, doppler, t_end_of_pass, signal_freq = pickle.load(f)
f.close()

print(len(doppler))

# add place for upsamp factor 2k = FS/FsDOp

#doppler = pchip_interpolate(range(0,(int(len(doppler)/FS_dop * FS)), 2000),doppler,range(0,(int(len(doppler)/FS_dop * FS)))) 
#doporig = doppler
#doppler = doppler[:int(len(doppler)/4)]
#doppler = 8000 - np.arange(0,len(doporig))*16000/len(doporig)  # this works

t1 = np.linspace(0,len(doppler), len(doppler))
t2 = np.linspace(0,len(doppler), int((len(doppler)*FS/FS_dop )/2))
print("Time", t1[-1], t2[-1])
#doppler = pchip_interpolate(t1,doppler,t2) 
ddop = int(FS/FS_dop)

doppler = np.concatenate((doppler,doppler[::-1]))
#doppler = pchip_interpolate(range(0,(int(len(doppler)/FS_dop * FS)), ddop),doppler,range(0,(int(len(doppler)/FS_dop * FS)))) 



# plt.figure(1)
# plt.plot(doppler)
# plt.xlabel('Doppler t = 1/100 s')
# plt.ylabel('Doppler')


# plt.figure(2)
# plt.plot(np.abs(np.diff(doppler())))
# plt.xlabel('Doppler t = 1/100 s')
# plt.ylabel('Diff(doppler)')
# plt.show()

print(len(doppler))
'''
midpoint = np.argwhere(np.diff(doppler)==min(np.diff(doppler)))[0]
doppler = doppler[int(midpoint-nums2):int(midpoint+nums2)]
doppler = np.concatenate((doppler,doppler[::-1]))
'''
'''
plt.figure(1)

plt.plot(doppler)
plt.plot(doporig)
plt.xlabel('Sample (1/FS)')
plt.ylabel('Frequency (Hz)')
plt.show()
'''
#t = np.linspace(0, len(doppler)/FS, len(doppler))
#nsp = 2**9*UPSAMP
#nsp = 2000
nsp = ddop
t = np.arange(0, nsp) / FS
# = np.linspace(0, (len(doppler)-1)/FS_dop, len(doppler))
#t_end = 20*60
#t = np.linspace(0, t_end, len(doppler))

'''
plt.figure(5)
plt.plot(np.diff(doppler))
plt.show()
'''
print(t[-1])
#doppler = 8000


output = []
phi_end=0
for dp in range(0,len(doppler)): 
    samp_shift= np.exp(1j * (2 * math.pi * doppler[dp] * t + phi_end))
    #samp_shift = samp_shift - phi_end
    #print(len(samp_shift))
    #print(len(samp_shift))
    output.extend(samp_shift)
    phi_end = np.angle(samp_shift[-1])
    print(phi_end)
#output = output.flatten()
print("OP len",len(output))
#print(len(doppler),len(t), len(output))

plt.figure(3)
plt.plot(np.abs(output[:3000]))
plt.show()
# plt.figure(2)
# plt.specgram(output)
# #plt.plot(np.angle(output))
# plt.show()
#wsoutput = output[
css_modulator = CssMod(N, SF, BW, FS, preamble, end_delimeter) 
bin_dat = np.float32(css_modulator.ang2bin_nopad(output))
bin_dat = bin_dat.tobytes()

file = open("GNURADIO_linear_8k_sep", 'bw')
file.write(bin_dat)
file.close()