'''
genExperimentFile :: Jackie Schellberg :: Last updated Feb 5, 2023 

For generating a single experiment, and keeping track of ground truth results. 
'''

import socket 
from css_mod import CssMod
from matplotlib import pyplot as plt
import numpy as np
import time
import math
import zmq
import random

#SF = 7 
SF = 9
N = 2**SF
UPSAMP = 10
NUMPKTS = 10
BW = 4

exp_root_folder = 'indoor_experiments'
exp_folder = '100pkts_sf9_bw'
fileout_name = "trial1"

symbols = []
symbols = np.ones(100)*5

preamble = [1,1]
end_delimeter = [3,3,3,3] 
css_modulator = CssMod(N, UPSAMP, preamble, end_delimeter) 
output = css_modulator.symbol2packet(symbols)

bin_dat = np.float32(css_modulator.ang2bin(output))

bsl = len(bin_dat)

bin_dat = np.append(np.float32(np.zeros((bsl))),bin_dat)
bin_dat = bin_dat.tobytes()

# input to gnuradio 
file = open(fileout_name, 'bw')
file.write(bin_dat)
file.close()

# save the ground truth symbols 
file = open(fileout_name + '_ref','w')
file.write(f"{symbols}")
file.close()
