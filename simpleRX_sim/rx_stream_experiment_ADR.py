'''
rx_stream.py :: Jackie Schellberg :: 10/19/2023 
Based off of code shared from Muhammad Osama Shahid. 

This is the top-level file to begin processing data obtained
from GNURADIO UDP SINK block. 
'''
import math
import multiprocessing
from multiprocessing import Manager
import socket

import time
from threading import Thread
from queue import Queue
import numpy as np
from utils import *
#from css_demod import CssDemod
from css_demod_2 import CssDemod
import pickle
import zmq

from matplotlib import pyplot as plt

#### Load data related to the experiment for BER/SNR/etc
#exp_root_folder = 'ber_desktop_testing'
# exp_root_folder = '../../experiment_data'
exp_root_folder = '../../experiment_data_synced'
exp_root_folder = '../../experiment_data_JUNE17'

f = open('../trial_under_test.txt')
exp_folder = f.read()
f.close()

# exp_folder = 'SF_7N_128BW_2500FS_200000NPKTS_5PLEN_100CR_0'

trial_name = "trial1"
#### 
# gnd_truth_data = '../simpleTX_sim/ground_truth'+'/'+exp_root_folder + '/' + exp_folder + '/' + trial_name + '.pkl'
# gnd_truth_data = exp_root_folder + '/' + exp_folder + '/' +"ground_truth_out.pkl"
# with open(gnd_truth_data,'rb') as f: 
    # SF, N, UPSAMP,NUMPKTS, BW, PAYLOAD_LEN,preamble,end_delimeter,symbols,CR = pickle.load(f)
SF = 7
N = 2 **7 
UPSAMP = 10
BW = 2500
PAYLOAD_LEN = 100
preamble = [0,0,0,0,0,0,0,0]
end_delimeter = [3,3,3,3] 
CR = 0
symbols = []
NUMPKTS = 0
############################### Variables ###############################
#RAW_FS = 450e3					# SDR's raw sampling freq
#RAW_FS = 200e3					# SDR's raw sampling freq
RAW_FS = 200000                # the queue size is selected so that no more than 1 packet may reside within a queue item
# RAW_FS = 3000000           # value should be kept <= expected length, so that we don't miss empty space
# RAW_FS=1250000
LORA_CHANNELS = [1]  # channels to process

UPSAMPLE_FACTOR = 4             		# factor to downsample to
#OVERLAP = 10 * UPSAMPLE_FACTOR#int(10 * RAW_FS / BW)
OVERLAP = 0
FS =200000
# FS =250000 # Min. BW for use with the rtlsdr
# Downsample the received signal to just the signal components of interest to lighten computation
GNURADIO_FS = 200000
DOPPLER_MAX = 18000
PRE_DOWNSAMP_FS =  (BW + DOPPLER_MAX)

PRE_DOWNSAMP = 1
# print(UPSAMP)
# UPSAMP = int(FS/BW)
# print("UPSAMP AFTER:", UPSAMP)
if (UPSAMP > 10):
    PRE_DOWNSAMP = int(FS/PRE_DOWNSAMP_FS) 
    # print(PRE_DOWNSAMP) 
    # some requirement for downsampling 
    # can't downsample by an invalid integer 
    while (FS % PRE_DOWNSAMP != 0):
        PRE_DOWNSAMP = PRE_DOWNSAMP - 1 
    # print("Initial Downsampling factor before processing: ", PRE_DOWNSAMP, "Into: ", FS/PRE_DOWNSAMP)
    # PRE_DOWNSAMP = 5
    FS = FS/PRE_DOWNSAMP 
else: 
    PRE_DOWNSAMP =1
    FS = FS/PRE_DOWNSAMP
# PRE_DOWNSAMP = 10
# UPSAMP = int(RAW_FS/BW)
UPSAMP = int(FS/BW)
# print("USAMP:", UPSAMP)
# UPSAMP = 10
PREAMBLE_SZ = 1*N*UPSAMP
# PREAMBLE_SZ = 1*N*(RAW_FS/BW)
# PREAMBLE_SZ = int(1*N*(RAW_FS/BW))
END_DELIMITER = end_delimeter
print("Downsampling the received signal: ", PRE_DOWNSAMP)
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
# FS = UPSAMP*BW
FS = UPSAMP*BW
# FS = RAW_FS
# FS = RAW_FS
# print("Sampling frequency:", FS)
#DB_THRESH = -7
# DB_THRESH = -11
# print("Bandwidth:",BW)
##########################################################################

def spawn_a_worker(my_channel, input_queue, output_queue):
    ###########################################################################################################
    # worker = YourDemodulatorImplementation(my_channel, input_queue, output_queue)
    # worker.start_consuming()
    ###########################################################################################################
    #css_demodulator = CssDemod(N, UPSAMP,PREAMBLE_SZ,END_DELIMITER,DB_THRESH);
    #css_demodulator = CssDemod(N, UPSAMP,PREAMBLE_SZ,END_DELIMITER,DB_THRESH, symbols,PAYLOAD_LEN,NUMPKTS,SF);
    css_demodulator = CssDemod(N, UPSAMP,PREAMBLE_SZ,END_DELIMITER, symbols,PAYLOAD_LEN,NUMPKTS,SF,BW,FS,CR);
    outfile = ('../simpleTX_sim/'+exp_root_folder + '/' + exp_folder + '/' + 'error_out')
    css_demodulator.setErrorMeasurementFile(outfile)
    print("Started")
    max_queue_cnt = 10
    while (True):
        if (input_queue.qsize() > 0): 
            [queue, queue_time] = input_queue.get()  
            # print("QUEUE_TIME:", queue_time)
            output = []
            css_demodulator.css_demod(my_channel, queue, output, queue_time)     
            if (len(output) >= 1):
                #print(output)
                print("====")
    print("Done.")
        
def IQ_SOURCE(chan, chan_num):
# def IQ_SOURCE(chan, chan_num, exp_time):
    context = zmq.Context() 
    socket = context.socket(zmq.SUB)
    socket.connect("tcp://127.0.0.1:55555")
    socket.setsockopt(zmq.SUBSCRIBE, b'')
    
    recver = list()
    counter = 0
    tmp = np.array([])
    atime = time.time()
    
    while True:
        if socket.poll(10) != 0:
            message = socket.recv() 
            
            if len(message) > 0:
                counter += len(message)
                data = np.frombuffer(message, dtype=np.complex64, count=-1) 
                recver.append(data)
                
                if counter > int(RAW_FS + OVERLAP * 2) * 8:  
                    final = np.concatenate(recver).ravel()
                    final = np.concatenate((tmp, final))                    # combine with previous
                    xl =len(tmp)
                    tmp = final[int(RAW_FS):]                               # get overlap
                    final = final[0:int(RAW_FS)]
                    
                    # decimate the signal to BW + DOPPLER MAX to quicken computation (especially in small BW case) 
                    # chan.put((final, time.time()))
                    chan.put((final[::PRE_DOWNSAMP], time.time()))
                    recver = list()
                    counter = len(tmp) * 8
                    atime += 1
                #break;
            else:
                time.sleep(0.1)

big_q = multiprocessing.Queue()
channel_streams = []
if __name__ == "__main__":
    # manager = Manager()
    
    for i in LORA_CHANNELS:
        in_queue = multiprocessing.Queue()
        channel_streams.append(in_queue)
        multiprocessing.Process(target=spawn_a_worker, args=(i, in_queue, big_q)).start()

    time.sleep(2.0)
    for i in range(len(LORA_CHANNELS)):
        #print(LORA_CHANNELS[i])
        # multiprocessing.Process(target=IQ_SOURCE, args=(channel_streams[i], LORA_CHANNELS[i])).start()
        multiprocessing.Process(target=IQ_SOURCE, args=(channel_streams[i], LORA_CHANNELS[i])).start()

    time.sleep(7260)