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
from pathlib import Path
import pickle 

#SF = 7 
SF = 9
N = 2**SF
UPSAMP = 10
NUMPKTS = 10
BW = 4
PAYLOAD_LEN = 100 
exp_root_folder = 'ber_desktop_testing'
exp_folder = '10pkts_sf9_20kbw_payload100'
trial_name = "trial1"

Path(exp_root_folder+'/'+exp_folder).mkdir(parents=True, exist_ok=True)
Path('ground_truth/'+exp_root_folder+'/'+exp_folder+'/').mkdir(parents=True, exist_ok=True)
fileout_name = exp_root_folder+'/'+exp_folder+'/'+trial_name
gt_fileout_name = 'ground_truth/'+exp_root_folder+'/'+exp_folder+'/'+trial_name+'.pkl'

symbols = []
#symbols = np.ones(100)*5
for i in range(0, PAYLOAD_LEN):
    symbols.append(random.randint(0,N))

preamble = [1,1]
end_delimeter = [3,3,3,3] 
css_modulator = CssMod(N, UPSAMP, preamble, end_delimeter) 
output = css_modulator.symbol2packet(symbols)

bin_dat = np.float32(css_modulator.ang2bin(output))

bsl = len(bin_dat)

#bin_dat = np.append(np.float32(np.zeros((bsl))),bin_dat)
#bin_dat = np.append(np.float32(np.zeros((bsl))),bin_dat)

bin_dat = np.tile(bin_dat,NUMPKTS)

# add zero padding so gnuradio has something to consume 
bin_dat = np.append(np.float32(np.zeros((bsl))),bin_dat)

bin_dat = bin_dat.tobytes()

# input to gnuradio 
file = open(fileout_name, 'bw')
file.write(bin_dat)
file.close()

# save the ground truth symbols to pkl 
with open(gt_fileout_name, 'wb') as f: 
    pickle.dump([SF, N, UPSAMP,NUMPKTS, BW, PAYLOAD_LEN,preamble,end_delimeter,symbols],f)
'''
file = open(fileout_name + '_ref','w')
file.write(f"{symbols}")
file.close()
'''