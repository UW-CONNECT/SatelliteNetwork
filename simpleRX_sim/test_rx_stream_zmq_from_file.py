import time
from threading import Thread
from queue import Queue
import numpy as np
from utils import *
from css_demod import CssDemod

import zmq

from matplotlib import pyplot as plt

# CSS Specifics to be known ahead of time
SF = 7 
N = 2**SF
UPSAMP = 10;
PREAMBLE_SZ = N*UPSAMP

END_DELIMITER = [3,3,3,3] 

# Threshold envelope; at what power level do we expect to see a packet arrive? 
# For low power scenario, this will have to be substituted 
#self.DB_THRESH = -13 # simulation with .005 noise voltage
#self.DB_THRESH = -25 # for simulation without any noise; should be perfect
#self.DB_THRESH = -27
#self.DB_THRESH = -20
#self.DB_THRESH = -7
#self.DB_THRESH = -30
DB_THRESH = -33.4
#self.DB_THRESH = -8   # sim at .035 noise


''' Old data using single packet size section ''' 
#filename =r'C:\Users\schellberg\Documents\schellberg\SatelliteNetwork\simpleRX_sim\rxAttenuator0dB'
#filename =r'C:\Users\schellberg\Documents\schellberg\SatelliteNetwork\simpleRX_sim\rxDataTmp\rxAttenuator0dB_2onespreamb_t5'
#filename =r'C:\Users\schellberg\Documents\schellberg\SatelliteNetwork\simpleRX_sim\rxDataTmp\rxAttenuator0dB_8onespreamb_t1'
#filename =r'C:\Users\schellberg\Documents\schellberg\SatelliteNetwork\simpleRX_sim\rxDataTmp\rxAttenuator0dB_8onespreamb_ALT_t1'

def load_file(file_path):
    fi_1 = open(file_path,'rb')
    x_inter = fi_1.read()
    x_inter_1=np.frombuffer(x_inter, dtype=np.float32)
    fi_1.close()
        
    x_1 = x_inter_1[::2] + 1j*x_inter_1[1::2] 
    return x_1
    
queue = load_file(filename)

#queue = CssDemod(N, UPSAMP,PREAMBLE_SZ,END_DELIMITER);
#queue = queue[11000:]
#queue = queue[1150000:]
print("Started")
max_queue_cnt = 10
#while (True):
css_demodulator = CssDemod(N, UPSAMP,PREAMBLE_SZ,END_DELIMITER);
output = []
css_demodulator.css_demod([], queue, output)     
if (len(output) > 1):
    print(output)
print("Done.")
    