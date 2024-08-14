import numpy as np 
import math
from matplotlib import pyplot as plt
def load_file(file_path):
    fi_1 = open(file_path,'rb')
    x_inter = fi_1.read()
    x_inter_1=np.frombuffer(x_inter, dtype=np.float32)
    fi_1.close()
        
    x_1 = x_inter_1[::2] + 1j*x_inter_1[1::2] 
    return x_1
FS = 200e3
  
doppler_out = load_file(r'C:\Users\Patron\Documents\zmqDoppler_2')
    
# doppler_out = load_file(r'J:\schellberg\indoor_exp_feb_2024\SatelliteNetwork-main\misc_test_scripts\50Hz_VCO_slope')
time_vec = np.arange(0, len(doppler_out))/FS

Freq = (np.angle(doppler_out[1:] * np.conjugate(doppler_out[0:-1])));
Freq = Freq / (2*math.pi) * (FS);
doppler_slope = np.diff(Freq) / (1/FS)
print("Slope:", np.mean(doppler_slope))
# plt.figure(1)
# plt.plot(time_vec[0:-1],Freq)
# plt.show()