''' Constant Frequency Slope :: Jackie Schellberg :: 6/24/2024 '''
import pickle 
from scipy.interpolate import pchip_interpolate
import numpy as np 
import math
from matplotlib import pyplot as plt
import sys
sys.path.insert(0, '../simpleTX_sim')
from css_mod import CssMod
# infile = "J:\schellberg\indoor_exp_feb_2024\doppler_simulation_files\ISS_doppler_sim_437_8.pkl"
# nums2 = FS # number of doppler samples (this corresponds to 2 seconds)
# filename = '50Hz_doppler_slope'
# info not needed
SF = 9 
preamble = [] 
end_delimeter = []
N = 2**SF
UPSAMP = 10
BW = 2500


FS =200000
SLOPE_TIME = 10*60*FS # number of samples, minutes * seconds * samples/second 

doppler_slope = 50 # Rate of change of doppler, Hz/second
# SLOPE_TIME = (FS/2)/doppler_slope * FS  # temporary 
# starting_frequency = 0
# starting_frequency = 0

time_vec = (np.arange(0, SLOPE_TIME) / FS )
# doppler = (np.arange(0, SLOPE_TIME) / FS ) * doppler_slope
doppler = (np.arange(0, SLOPE_TIME) / FS ) * doppler_slope
# doppler = 

# plt.figure(1)
# plt.plot(time_vec, doppler)
# plt.title("Doppler")
# plt.show()


filename = '../../../'+ str(doppler_slope) + 'Hz_doppler_slope'

# plt.figure(1)
# plt.plot(doppler)
# plt.show()
# output = np.exp(1j*( 2*math.pi*(doppler*time_vec) ))
output = np.exp(1j*(2*math.pi * doppler_slope * (time_vec**2) / 2))

Freq = (np.angle(output[1:] * np.conjugate(output[0:-1])));
Freq = Freq / (2*math.pi) * (FS);
# plt.figure(1)
# plt.plot(Freq)    # Sanity check that doppler vector is producing the expected outputs.
# plt.show()
# for dp in range(1,len(doppler)): 
    # # samp_shift= np.exp(1j * (2 * math.pi * doppler[dp] * t + phi_end))
    # k = (doppler[dp] - doppler[dp-1])/len(t)
    # samp_shift = np.exp(1j*(phi_end + 2*math.pi*(t*(doppler[dp-1] +0.5*k*t)) ))
    # #samp_shift = samp_shift - phi_end
    # #print(len(samp_shift))
    # #print(len(samp_shift))
    # output.extend(samp_shift)
    # doppler_ref.append(np.tile(doppler[dp], len(t)))
    # phi_end = np.angle(samp_shift[-1])
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

# with open(filename + '_time.pkl', 'wb') as f:  # Python 3: open(..., 'wb')
    # pickle.dump([doppler_ref, FS], f)

# plt.figure(1)
# plt.plot(doppler_ref) 
# plt.show()    
    
# plt.figure(3)
# plt.plot(np.abs(output[:3000]))
# plt.show()
# plt.figure(2)
# plt.specgram(output)
# plt.plot(np.angle(output))
# plt.show()
#wsoutput = output[
CR = 0
css_modulator = CssMod(N, SF, BW, FS, preamble, end_delimeter, CR) 
bin_dat = np.float32(css_modulator.ang2bin_nopad(output))
bin_dat = bin_dat.tobytes()

file = open(filename, 'bw')
file.write(bin_dat)
file.close()
