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
infile = "J:\schellberg\indoor_exp_feb_2024\doppler_simulation_files\ISS_doppler_sim_437_8.pkl"
nums2 = FS # number of doppler samples (this corresponds to 2 seconds)
filename = 'doppler_sync_testing'
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

# bin_dat = np.append(bin_dat,np.float32(np.zeros((buffer_dat*8)))) # for simulation in GNURADIO
            
# nsamps_total = len(bin_dat) / (2)
# bin_dat = np.float32(doppler)
# bin_dat = bin_dat.tobytes()
# f = open("doppler_out", "wb")
# f.write(bin_dat)
# f.close()

# plt.figure(1)
# plt.plot(doppler)
# plt.show()

# mid = int(len(doppler)/2) - 300
mid = 42300
nsamps = 100 * 60 * 1
# doppler = doppler[mid-nsamps:mid+nsamps]
# doppler = np.ones(8000)*8000

# add place for upsamp factor 2k = FS/FsDOp

#doppler = pchip_interpolate(range(0,(int(len(doppler)/FS_dop * FS)), 2000),doppler,range(0,(int(len(doppler)/FS_dop * FS)))) 
# doporig = doppler
# doppler = doppler[:int(len(doppler)/4)]
# doporig = np.ones(3000)
# doppler = 8000 - np.arange(0,len(doporig))*16000/len(doporig)  # this works

# t1 = np.linspace(0,len(doppler), len(doppler))
# t2 = np.linspace(0,len(doppler), int((len(doppler)*FS/FS_dop )/2))
# print("Time", t1[-1], t2[-1])
# #doppler = pchip_interpolate(t1,doppler,t2) 
ddop = int(FS/FS_dop)

dop_ax = np.arange(0, len(doppler)) / (FS_dop)

sampling_interval = 1
dop_ax = dop_ax[::sampling_interval]
doppler = doppler[::sampling_interval]

dt = sampling_interval / FS_dop
print("DT: ", dt)

doppler_overhead = dop_ax[np.argmin(np.abs(doppler))]
# doppler_overhead = np.argmin(np.abs(doppler))
# print(doppler_overhead)
plt.show()
plt.figure(1)
plt.plot(dop_ax, doppler)
# plt.plot( doppler)
plt.xlabel('Time (s)')
plt.ylabel('Doppler Shift (Hz)')
plt.axvline(doppler_overhead,color='r')
plt.title('ISS Pass Doppler Shift ')

plt.figure(2)
plt.plot(dop_ax[0:-1], np.abs(np.diff(doppler)/dt))
# augmented_dop = np.concatenate(( doppler, [doppler[-1]]))
# plt.plot( np.abs(np.diff(augmented_dop)))
plt.xlabel('Time (s)')
plt.ylabel('Abs. Value of Doppler Shift Change (Hz)')
plt.axvline(doppler_overhead,color='r')
plt.title('ISS Pass Doppler Shift Rate of Change (Slope)')

plt.figure(3) 
plt.plot(dop_ax[1:-1], np.abs(np.diff(doppler,2) /dt**2))
# plt.plot( np.abs(np.diff(doppler,2 )))
# augmented_dop = np.concatenate(( doppler, [doppler[-1]], [doppler[-1]]))

# plt.figure(6)
# plt.plot(augmented_dop)
# plt.show()

# plt.plot( np.abs(np.diff(augmented_dop,2 )))
plt.xlabel('Time (s)')
plt.ylabel('Abs. Value of Doppler Slope Change (Hz)')
plt.axvline(doppler_overhead,color='r')
plt.title('ISS Pass Doppler Shift Acceleration')


plt.show()



# doppler = np.concatenate((doppler,doppler[::-1]))
#doppler = pchip_interpolate(range(0,(int(len(doppler)/FS_dop * FS)), ddop),doppler,range(0,(int(len(doppler)/FS_dop * FS)))) 



# plt.figure(1)
# plt.plot(doppler)
# plt.xlabel('Doppler t = 1/100 s')
# plt.ylabel('Doppler')
# plt.show()

# plt.figure(2)
# plt.plot(np.abs(np.diff(doppler())))
# plt.xlabel('Doppler t = 1/100 s')
# plt.ylabel('Diff(doppler)')
# plt.show()

# print(len(doppler))
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
print("Num samples:", nsp)
t = np.arange(0, nsp) / FS
# = np.linspace(0, (len(doppler)-1)/FS_dop, len(doppler))
#t_end = 20*60
#t = np.linspace(0, t_end, len(doppler))

'''
plt.figure(5)
plt.plot(np.diff(doppler))
plt.show()
'''
# print(t[-1])
#doppler = 8000


output = []
phi_end=0
doppler_ref = [] # keep the upsampled doppler vector 
# for dp in range(0,len(doppler)): 
for dp in range(1,len(doppler)): 
    # samp_shift= np.exp(1j * (2 * math.pi * doppler[dp] * t + phi_end))
    k = (doppler[dp] - doppler[dp-1])/len(t)
    samp_shift = np.exp(1j*(phi_end + 2*math.pi*(t*(doppler[dp-1] +0.5*k*t)) ))
    #samp_shift = samp_shift - phi_end
    #print(len(samp_shift))
    #print(len(samp_shift))
    output.extend(samp_shift)
    doppler_ref.append(np.tile(doppler[dp], len(t)))
    phi_end = np.angle(samp_shift[-1])
    # print(phi_end)
#output = output.flatten()
print("OP len",len(output))
#print(len(doppler),len(t), len(output))

# plt.figure(1)
# plt.plot(output[:len(t)*2])
# plt.show()
output_time_vec = len(output)
# file = open(filename + '_time', 'bw')
# file.write(output_time_vec)
# file.close()

with open(filename + '_time.pkl', 'wb') as f:  # Python 3: open(..., 'wb')
    pickle.dump([doppler_ref, FS], f)

plt.figure(1)
plt.plot(doppler_ref) 
plt.show()    
    
# plt.figure(3)
# plt.plot(np.abs(output[:3000]))
# plt.show()
# plt.figure(2)
# plt.specgram(output)
# # plt.plot(np.angle(output))
# plt.show()
#wsoutput = output[
CR = 0
css_modulator = CssMod(N, SF, BW, FS, preamble, end_delimeter, CR) 
bin_dat = np.float32(css_modulator.ang2bin_nopad(output))
bin_dat = bin_dat.tobytes()

file = open(filename, 'bw')
file.write(bin_dat)
file.close()