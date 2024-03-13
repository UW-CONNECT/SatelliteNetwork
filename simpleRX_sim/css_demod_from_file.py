import time
from threading import Thread
from queue import Queue
import numpy as np
from utils import *
from css_demod_2 import CssDemod
import zmq

from matplotlib import pyplot as plt

# CSS Specifics to be known ahead of time
SF = 7
#SF = 10
N = 2**SF
UPSAMP = 10;
PREAMBLE_SZ = 1*N*UPSAMP
#PREAMBLE_SZ = 3*N*UPSAMP

END_DELIMITER = [3,3,3,3] 

# Threshold envelope; at what power level do we expect to see a packet arrive? 
# For low power scenario, this will have to be substituted 
DB_THRESH = -8   # sim at .035 noise

#filename =r'C:\Users\schellberg\Documents\schellberg\SatelliteNetwork\simpleTX_sim\test_files\SF9100s_035chan'
#filename = r'C:\Users\schellberg\Documents\schellberg\Standard_LoRa\test_3ones'
# filename = r'J:\schellberg\indoor_exp_feb_2024\SatelliteNetwork-main\simpleTX_sim\two_usrp_testing\sf9_20k_preamb3_ones_3\test_3ones_2_LOWNOISE'
#filename =  r'J:\schellberg\indoor_exp_feb_2024\experiment_data\SF_7N_TX'
filename = r'J:\schellberg\indoor_exp_feb_2024\simulated_data_SNR\sf7_noDop_SNR1-5hundredth_4dB'
#filename = r'J:\schellberg\spurious_906MHz_10pkts_SF9' #incorrect preamble here
#filename = r'J:\schellberg\indoor_exp_feb_2024\SatelliteNetwork-main\simpleTX_sim\two_usrp_testing\sf9_20k_preamb3_ones_3\test_3ones_2_HIGHSNR'
def load_file(file_path):
    fi_1 = open(file_path,'rb')
    x_inter = fi_1.read()
    x_inter_1=np.frombuffer(x_inter, dtype=np.float32)
    fi_1.close()
        
    x_1 = x_inter_1[::2] + 1j*x_inter_1[1::2] 
    return x_1
    
queue = load_file(filename)
# queue = queue[:int(len(queue)/9) * 2] #4e5 expected start idx
print(len(queue))
plt.figure(1)
plt.plot(queue)
plt.show()
print("Started")
max_queue_cnt = 10
#while (True):
symbols= np.ones(100)*5
PAYLOAD_LEN= 100
NUMPKTS=1
BW = 20000
FS = 200000
css_demodulator = CssDemod(N, UPSAMP,PREAMBLE_SZ,END_DELIMITER,DB_THRESH, symbols,PAYLOAD_LEN,NUMPKTS,SF,BW,FS);
output = []
css_demodulator.css_demod([], queue, output)     

print(output)
print("Done.")
    