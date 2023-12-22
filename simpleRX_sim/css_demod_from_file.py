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
DB_THRESH = -8   # sim at .035 noise

filename =r'C:\Users\schellberg\Documents\schellberg\SatelliteNetwork\simpleTX_sim\thousand5_tests\SF9_035chan'

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
css_demodulator = CssDemod(N, UPSAMP,PREAMBLE_SZ,END_DELIMITER,DB_THRESH);
output = []
css_demodulator.css_demod([], queue, output)     
if (len(output) > 1):
    print(output)
print("Done.")
    