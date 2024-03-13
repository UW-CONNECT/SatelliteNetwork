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


SF = 12
#SF = 9
N = 2**SF
UPSAMP = 10
NUMPKTS = 1
BW = 20000 
FS = 200000
PAYLOAD_LEN = 100
exp_root_folder = '../../experiment_data'
#exp_root_folder = '../../simulated_data'
#exp_folder = str(NUMPKTS) + '_pkts_' + str(SF) + '_SF'  
exp_folder = 'SF_' + str(SF) + 'N_' + str(N) + 'BW_' + str(BW) + 'FS_' + str(FS) +\
'NPKTS_' +str(NUMPKTS) + 'PLEN_' +str(PAYLOAD_LEN)
trial_name = "trial1"

Path(exp_root_folder+'/'+exp_folder).mkdir(parents=True, exist_ok=True)
Path('ground_truth/'+exp_root_folder+'/'+exp_folder+'/').mkdir(parents=True, exist_ok=True)
fileout_name = exp_root_folder+'/'+exp_folder+'/'+trial_name
gt_fileout_name = exp_root_folder+'/'+exp_folder+'/'+"ground_truth_out.pkl"
#gt_fileout_name = 'ground_truth/'+exp_root_folder+'/'+exp_folder+'/'+trial_name+'.pkl'
#symbols = []
#symbols = np.ones(PAYLOAD_LEN)*5

#symbols = []
#symbols = np.ones(PAYLOAD_LEN)*5

symbols = []
#for ss in range(0,1000):
#    symbols.append(random.randint(1,N-1))
for ss in range(0,100):
    symbols.append(random.randint(1,N-1))
#symbols = np.ones(100)*5
#symbols = [1]

#preamble = [0,0,0,0,0,0,0,0]
#end_delimeter = [0,0,0,0,0,0,0,0] 
#preamble = [1, 0,0,1,1,0,0,1]
#preamble = [1,1]
preamble = [1,1,1]
#preamble = [1,1,1,1]
#end_delimeter = [1,0,5,0,0,10,1,1] 
end_delimeter = [3,3,3,3] 
css_modulator = CssMod(N, SF, BW, FS, preamble, end_delimeter) 
output = css_modulator.symbol2packet(symbols)

print("OUTPUT len", len(output))

bin_dat = np.float32(css_modulator.ang2bin(output))
#print(bin_dat)


bsl = len(bin_dat)
#bin_dat = np.concatenate((bin_dat,bin_dat,bin_dat,bin_dat,bin_dat))

#print(bin_dat.shape)
#bin_dat = np.concatenate((np.zeros((bsl)), bin_dat))
#print(bin_dat.shape)
#print(len(bin_dat))
buffer_dat = N*UPSAMP*5*32
bin_dat = np.append(np.float32(np.zeros((buffer_dat))),bin_dat)

bin_dat = np.tile(bin_dat, NUMPKTS)
#bin_dat = np.append(np.float32(np.zeros((bsl))),bin_dat)
bin_dat = bin_dat.tobytes()
'''
symbols = []
for i in range(0, PAYLOAD_LEN):
    symbols.append(random.randint(0,N-1))
#symbols = np.asarray(symbols)
#print(len(symbols),len(symbols2))
preamble = [1,1,1]
end_delimeter = [3,3,3,3] 
css_modulator = CssMod(N, UPSAMP, preamble, end_delimeter) 
output = css_modulator.symbol2packet(symbols)

bin_dat = np.float32(css_modulator.ang2bin(output))

bsl = len(bin_dat)

#bin_dat = np.append(np.float32(np.zeros((bsl))),bin_dat)
#bin_dat = np.append(np.float32(np.zeros((bsl))),bin_dat)

# add zero padding so gnuradio has something to consume 
buffer_dat = N*UPSAMP*5*32
bin_dat = np.append(np.float32(np.zeros((buffer_dat))),bin_dat)

#bin_dat = np.append(np.float32(np.zeros((bsl))),bin_dat)
bin_dat = bin_dat.tobytes()
#bin_dat = np.tile(bin_dat,NUMPKTS)
'''



# input to gnuradio 
print(fileout_name)
print(len(bin_dat))
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