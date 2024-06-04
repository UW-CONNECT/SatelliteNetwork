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


# SF = 7
# SF_list = [7,8,9]
SF_list = [7, 8, 9]
#SF = 9

NUMPKTS = 25
# BW = 20000 
# BW_list = [5000, 10000, 20000, 50000]
# BW_list = [5000, 10000, 20000, 40000, 50000]
# BW_list = [2500, 10000, 20000]
BW_list = [2500, 20000]
# BW_list = [3000,5000,10000,20000]
# BW_list = [20000]
# print(len(BW_list))
# BW = 125000 
# FS = 200000
FS = 200000
# FS = 1000000
# FS = BW*10
PAYLOAD_LEN = 100
# CR_LIST = [0,1,2,3,4] 
# CR_LIST = [0,3] 
CR_LIST = [0]

symbols = []
#for ss in range(0,1000):
#    symbols.append(random.randint(1,N-1))

# preamble = [0,0,0,0,0,0,0,0,0,0]
preamble = [0,0,0,0,0,0,0,0]
end_delimeter = [3,3,3,3] 

trial_num_name_out = 0
for SF in SF_list:
    N = 2**SF
    symbols = []
    for ss in range(0,PAYLOAD_LEN):
        symbols.append(random.randint(1,N-1))

    for BW in BW_list:
        for CR in CR_LIST:
            
            UPSAMP = 10 # grandfathered
        
        # CR = 0 # coding rate, 0 for no hamming encoding, 1-4 for other corresponding options
            # exp_root_folder = '../../experiment_data'
            exp_root_folder = '../../experiment_data_june'
            #exp_root_folder = '../../simulated_data'
            #exp_folder = str(NUMPKTS) + '_pkts_' + str(SF) + '_SF'  
            exp_folder = 'SF_' + str(SF) + 'N_' + str(N) + 'BW_' + str(BW) + 'FS_' + str(FS) +\
            'NPKTS_' +str(NUMPKTS) + 'PLEN_' +str(PAYLOAD_LEN) +'CR_' +str(CR)
            # exp_folder = str(trial_num_name_out)
            trial_name = "trial1"

            Path(exp_root_folder+'/'+exp_folder).mkdir(parents=True, exist_ok=True)
            Path('ground_truth/'+exp_root_folder+'/'+exp_folder+'/').mkdir(parents=True, exist_ok=True)
            fileout_name = exp_root_folder+'/'+exp_folder+'/'+trial_name
            # fileout_name = exp_root_folder+'/'+str(trial_num_name_out)+'/'+trial_name
            # fileout_name = exp_root_folder+'/'+exp_folder+'/'+
            # trial_num_name_out = trial_num_name_out + 1
            trial_num_name_out = trial_num_name_out + 1
            gt_fileout_name = exp_root_folder+'/'+exp_folder+'/'+"ground_truth_out.pkl"
            
            css_modulator = CssMod(N, SF, BW, FS, preamble, end_delimeter, CR) 
            output = css_modulator.symbol2packet(symbols)

            print("OUTPUT len", len(output))

            bin_dat = np.float32(css_modulator.ang2bin(output))
            #print(bin_dat)


            bsl = len(bin_dat)
            buffer_dat = N*UPSAMP*5*32
            bin_dat = np.append(np.float32(np.zeros((buffer_dat))),bin_dat)

            bin_dat = np.tile(bin_dat, NUMPKTS)
            bin_dat = np.append(bin_dat,np.float32(np.zeros((buffer_dat*8)))) # for simulation in GNURADIO
            
            nsamps_total = len(bin_dat) / (2)
            
            bin_dat = bin_dat.tobytes()
            
            
            DATA_FILE_PROC_TIME = nsamps_total/FS
            data_proc_time_filename = exp_root_folder+'/'+exp_folder+'/'+"proc_file_time"
            f = open(data_proc_time_filename, "w")
            f.write(str(DATA_FILE_PROC_TIME))
            f.close()
            
            
            
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
                pickle.dump([SF, N, UPSAMP,NUMPKTS, BW, PAYLOAD_LEN,preamble,end_delimeter,symbols, CR],f)
            '''
            file = open(fileout_name + '_ref','w')
            file.write(f"{symbols}")
            file.close()
            '''