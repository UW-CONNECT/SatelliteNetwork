'''

Checks doppler file timestamp and determines current doppler shift. Assigns 
SF based on the rate-of-change of doppler. 

[TODO] What to do if we have a TLE rather than a doppler start time? 
-- Replace this code with something to check 

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

# GNURADIO And signal params 
FS = 200000
BW = 2500
preamble = [0,0,0,0,0,0,0,0]
end_delimeter = [3,3,3,3] 
CR = 0 # no error correction for now 
PAYLOAD_LEN = 100 # number of symbols within a payload
exp_root_folder = '../../experiment_data_ADR'
exp_folder = 0 # counts upwards with every packet sent 
NUMPKTS = 1
# TODO; save the reference doppler and the time array separately, and interpolate these 
doppler_time_file = r'J:\schellberg\indoor_exp_feb_2024\SatelliteNetwork-main\simpleTX_sim\doppler_sim_scripts\doppler_sync_testing_time.pkl'
live_doppler_file = r'J:\schellberg\indoor_exp_feb_2024\SatelliteNetwork-main\transponder_sim\doppler_start_time_jun17'

# output_time_vec = np.arange(0, len(output))/ FS 
f = open(doppler_time_file,'rb')
[doppler_ref, FS_dop] = pickle.load(f)          # [TODO] why does this take so long to load? 

doppler_ref = np.array(doppler_ref)
doppler_ref = doppler_ref[:,0]
f.close()

f = open(live_doppler_file, 'r')
transponder_start_time = float(f.read())
f.close()

print("Linear transponder with doppler start time: ", transponder_start_time)

# print(len(FS_dop)) 
FS_dop = 100
doppler_time_vec = np.arange(0, len(doppler_ref))/FS_dop
# print(doppler_time_vec[-1],len(doppler_ref),len(doppler_ref)/FS_dop)

# plt.figure(1)
# plt.plot(doppler_time_vec, doppler_ref)
# plt.show()

doppler_1s = doppler_ref[::FS_dop]
diff_dop_Hz = np.abs(np.diff(doppler_1s))
# plt.hist(diff_dop_Hz) 
# plt.xlabel('Doppler Change (Hz/s)')
# plt.ylabel('Duration/Frequency')
# plt.show()
# Consider the time of the longest chirp; this is the rate at which we will transmit packets 
TX_TIME = 2

# Check the rate-of-change of doppler and map to an SF [TODO]  
DDOP_SF = {}

# Read the doppler start time file [TODO]
# doppler_start_time = 0 

def get_rate(d_doppler):
    '''
    Uses the rate of change of doppler, Hz/s, to determine the rate SF={7,8,9,10,11,12}
    '''
    SF = 12 # default most reliability 
    
    if 80 <= d_doppler  :
        SF = 7
    elif 60 <= d_doppler < 80 :
        SF = 8
    elif 40 <= d_doppler < 60 :
        SF = 9
    elif 20 <= d_doppler < 40 :
        SF = 10
    elif 10 <= d_doppler < 20 :
        SF = 11 
    elif d_doppler < 10 :
        SF = 12 
    else :
        print("Invalid doppler rate. ")
        return 0 
    return SF

def gen_pkt(exp_root_folder,exp_folder, SF,BW,FS,preamble, end_delimeter,CR):
    '''
    Generate a randomized packet at the specified SF
    '''            
    UPSAMP = int(FS/BW)
    N = 2**SF
    symbols = []
    for ss in range(0,PAYLOAD_LEN):
        symbols.append(random.randint(1,N-1))
    # UPSAMP = 10 # grandfathered

    Path(exp_root_folder+'/'+exp_folder).mkdir(parents=True, exist_ok=True)
    Path('ground_truth/'+exp_root_folder+'/'+exp_folder+'/').mkdir(parents=True, exist_ok=True)
    gt_fileout_name = exp_root_folder+'/'+exp_folder+'/'+"ground_truth_out.pkl"
    
    css_modulator = CssMod(N, SF, BW, FS, preamble, end_delimeter, CR) 
    output = css_modulator.symbol2packet(symbols)

    # print("OUTPUT len", len(output))

    bin_dat = np.float32(css_modulator.ang2bin(output))
    #print(bin_dat)


    bsl = len(bin_dat)
    buffer_dat = N*UPSAMP*5*32
    bin_dat = np.append(np.float32(np.zeros((buffer_dat))),bin_dat)

    bin_dat = np.tile(bin_dat, NUMPKTS)
    bin_dat = np.append(bin_dat,np.float32(np.zeros((buffer_dat*8)))) # for simulation in GNURADIO
    
    nsamps_total = len(bin_dat) / (2)
    
    bin_dat = bin_dat.tobytes()
    
    # save the ground truth symbols to pkl 
    with open(gt_fileout_name, 'wb') as f: 
        pickle.dump([SF, N, UPSAMP,NUMPKTS, BW, PAYLOAD_LEN,preamble,end_delimeter,symbols, CR],f)

    return bin_dat
# def send_samples_ZMQ(data):
    # '''
    # Send I/Q samples to ZMQ for use in GNURADIO 
    # '''
    
    # # experiment_start_time = time.time()
    
    # # experiment_time_path= EXP_ROOT_FOLDER + '\\' + str(TRIAL_FOLDER) + '\\' + "posix_start_time"
    # # f = open(experiment_time_path, "w")
    # # f.write(str(experiment_start_time))
    # # f.close()
    
    # # print("Started sending TX data over TCP. Start time: ", experiment_start_time)
    # # with open(FILE_PATH, "rb") as f:
        # # nbytes= 200000
        # # print(nbytes)
        # # while (byte := f.read(nbytes)):
            # # Do stuff with byte.        
    # zsocket.send(data)
            # # byte = np.frombuffer(byte, dtype=np.float32)
            # # print("Sent byte")
            # # time.sleep(.1)
    # print("Done sending the data over ZMQ")
context = zmq.Context() 
zsocket = context.socket(zmq.PUB) 
zsocket.bind("tcp://127.0.0.1:4444")    
print(doppler_ref.shape)   
doppler_prev = doppler_ref[-1]    
while (True):
    current_time = (time.time() - transponder_start_time) % doppler_time_vec[-1]   
    
    current_doppler = np.argmin(np.abs(np.subtract(doppler_time_vec, current_time)))
    # current_doppler = current_doppler[0]
    # print(current_doppler, sop)
    print("Live transponder doppler value: ", doppler_ref[current_doppler], "at time: ", current_time)
    current_doppler_value = doppler_ref[current_doppler]
    ddop = np.abs(current_doppler_value - doppler_prev) / TX_TIME 
    doppler_prev = current_doppler_value
    print("Current doppler rate: ", ddop,"Hz/s")
    print(  "Selected SF: ", get_rate(ddop))
    SF = get_rate(ddop)
    ''' Transmit the packet via ZMQ'''
    zsocket.send(gen_pkt(exp_root_folder,str(exp_folder), SF,BW,FS,preamble, end_delimeter,CR))
    exp_folder = exp_folder + 1
    
    # wait to Tx next packet 
    time.sleep(TX_TIME) # [TODO] Change this to match the period of the transmitted packet, so that we can Tx faster