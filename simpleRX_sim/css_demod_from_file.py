import time
from threading import Thread
from queue import Queue
import numpy as np
from utils import *
from css_demod_2 import CssDemod
import zmq

from matplotlib import pyplot as plt
import math
import multiprocessing
from multiprocessing import Manager
import socket

import time

import pickle

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
filename = r'C:\Users\schellberg\Documents\schellberg\Standard_LoRa\sf7_with_doppler'
#filename = r'J:\schellberg\spurious_906MHz_10pkts_SF9' #incorrect preamble here
#filename = r'J:\schellberg\indoor_exp_feb_2024\SatelliteNetwork-main\simpleTX_sim\two_usrp_testing\sf9_20k_preamb3_ones_3\test_3ones_2_HIGHSNR'


#### Load data related to the experiment for BER/SNR/etc
#exp_root_folder = 'ber_desktop_testing'
exp_root_folder = '../../experiment_data'

exp_folder = 'SF_7N_128BW_40000FS_200000NPKTS_100PLEN_100CR_0'

trial_name = "trial1"
#### 
gnd_truth_data = '../simpleTX_sim/ground_truth'+'/'+exp_root_folder + '/' + exp_folder + '/' + trial_name + '.pkl'
gnd_truth_data = exp_root_folder + '/' + exp_folder + '/' +"ground_truth_out.pkl"
with open(gnd_truth_data,'rb') as f: 
    SF, N, UPSAMP,NUMPKTS, BW, PAYLOAD_LEN,preamble,end_delimeter,symbols,CR = pickle.load(f)

def load_file(file_path):
    fi_1 = open(file_path,'rb')
    x_inter = fi_1.read()
    x_inter_1=np.frombuffer(x_inter, dtype=np.float32)
    fi_1.close()
        
    x_1 = x_inter_1[::2] + 1j*x_inter_1[1::2] 
    return x_1
    
queue = load_file(filename)
# queue = queue[:int(len(queue)/9) * 2] #4e5 expected start idx
# print(len(queue))
# plt.figure(1)
# plt.plot(queue)
# plt.show()
# print("Started")
# max_queue_cnt = 10
# while (True):
# symbols= np.ones(100)*5
# PAYLOAD_LEN= 100
# NUMPKTS=1
# BW = 20000
# FS = 200000
RAW_FS = 1000000                # the queue size is selected so that no more than 1 packet may reside within a queue item
# RAW_FS = 1000000           # value should be kept <= expected length, so that we don't miss empty space
# RAW_FS=1250000
LORA_CHANNELS = [1]  # channels to process

UPSAMPLE_FACTOR = 4             		# factor to downsample to
#OVERLAP = 10 * UPSAMPLE_FACTOR#int(10 * RAW_FS / BW)
OVERLAP = 0

# CSS Specifics to be known ahead of time
#SF = 9
#N = 2**SF
#UPSAMP = 10;
#print(SF, N, UPSAMP)
#PREAMBLE_SZ = int(len(preamble)/2)*N*UPSAMP
#PREAMBLE_SZ = 3*N*UPSAMP
FS = 200000
# UPSAMP = int(RAW_FS/BW)
UPSAMP = int(FS/BW)
PREAMBLE_SZ = 1*N*UPSAMP
# PREAMBLE_SZ = 1*N*(RAW_FS/BW)
# PREAMBLE_SZ = int(1*N*(RAW_FS/BW))
END_DELIMITER = end_delimeter

# Threshold envelope; at what power level do we expect to see a packet arrive? 
# For low power scenario, this will have to be substituted 
#self.DB_THRESH = -13 # simulation with .005 noise voltage
#DB_THRESH = -25 # for simulation without any noise; should be perfect
#self.DB_THRESH = -27
#self.DB_THRESH = -20
#self.DB_THRESH = -7
#self.DB_THRESH = -30
#DB_THRESH = -33.4
#DB_THRESH = -8   # sim at .035 noise
#DB_THRESH = -5.5
FS = UPSAMP*BW
# print("Sampling frequency:", FS)
#DB_THRESH = -7
DB_THRESH = -11

css_demodulator = CssDemod(N, UPSAMP,PREAMBLE_SZ,END_DELIMITER,DB_THRESH, symbols,PAYLOAD_LEN,NUMPKTS,SF,BW,FS,CR);
output = []
css_demodulator.css_demod([], queue, output)     

print(output)
print("Done.")
    