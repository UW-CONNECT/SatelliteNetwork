import time
from threading import Thread
from queue import Queue
import numpy as np
from utils import *
from css_demod_2 import CssDemod
# css_demod import CssDemod
import zmq
import pickle
from matplotlib import pyplot as plt

# CSS Specifics to be known ahead of time
SF = 9
N = 2**SF
UPSAMP = 10;
PREAMBLE_SZ = N*UPSAMP

'''
What is the preamble size? 
'''

END_DELIMITER = [3,3,3,3] 

# Threshold envelope; at what power level do we expect to see a packet arrive? 
# For low power scenario, this will have to be substituted 
DB_THRESH = -8   # sim at .035 noise
DB_THRESH = -9.5
#filename =r'C:\Users\schellberg\Documents\schellberg\SatelliteNetwork\simpleTX_sim\ber_desktop_testing\'
filename = r'C:\Users\schellberg\Documents\schellberg\SatelliteNetwork\simpleTX_sim\ber_desktop_testing\10pkts_sf9_20kbw_payload100\trial1_035chan'
def load_file(file_path):
    fi_1 = open(file_path,'rb')
    x_inter = fi_1.read()
    x_inter_1=np.frombuffer(x_inter, dtype=np.float32)
    fi_1.close()
        
    x_1 = x_inter_1[::2] + 1j*x_inter_1[1::2] 
    return x_1
    
queue = load_file(filename)

print("Started")
max_queue_cnt = 10
#while (True):
symbols= np.ones(100)*5
PAYLOAD_LEN= 100
NUMPKTS=1


exp_root_folder = 'ber_desktop_testing'
exp_folder = '10pkts_sf9_20kbw_payload100'
#exp_folder = '10000pkts_sf9_20kbw_payload100'
trial_name = "trial1"
gnd_truth_data = '../simpleTX_sim/ground_truth'+'/'+exp_root_folder + '/' + exp_folder + '/' + trial_name + '.pkl'
with open(gnd_truth_data,'rb') as f: 
    SF, N, UPSAMP,NUMPKTS, BW, PAYLOAD_LEN,preamble,end_delimeter,symbols = pickle.load(f)


css_demodulator = CssDemod(N, UPSAMP,PREAMBLE_SZ,END_DELIMITER,DB_THRESH, symbols,PAYLOAD_LEN,NUMPKTS);
output = []
#queue = queue[:1850000]
#queue = queue[:1850000*2]
css_demodulator.css_demod([], queue, output)     

print(output)
print("Done.")
    
    