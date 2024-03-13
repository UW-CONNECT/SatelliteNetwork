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
exp_root_folder = '../../experiment_data'
#exp_root_folder = '../../simulated_data'
exp_folder = 'SF_7N_128BW_20000FS_200000NPKTS_10PLEN_100'
#exp_folder = 'SF_8N_256BW_20000FS_200000NPKTS_10PLEN_100'
#exp_folder = 'SF_9N_512BW_20000FS_200000NPKTS_10PLEN_100'
#exp_folder = 'SF_10N_1024BW_20000FS_200000NPKTS_10PLEN_100'
#exp_folder = 'SF_9N_512BW_20000FS_200000NPKTS_100PLEN_100'
#exp_folder = 'SF_11N_2048BW_20000FS_200000NPKTS_10PLEN_100' 
#exp_folder = 'SF_12N_4096BW_20000FS_200000NPKTS_10PLEN_100'
#exp_folder = '10000pkts_sf9_20kbw_payload100'
trial_name = "trial1"
#### 
gnd_truth_data = '../simpleTX_sim/ground_truth'+'/'+exp_root_folder + '/' + exp_folder + '/' + trial_name + '.pkl'
gnd_truth_data = exp_root_folder + '/' + exp_folder + '/' +"ground_truth_out.pkl"
with open(gnd_truth_data,'rb') as f: 
    SF, N, UPSAMP,NUMPKTS, BW, PAYLOAD_LEN,preamble,end_delimeter,symbols = pickle.load(f)
############################### Variables ###############################
#RAW_FS = 450e3					# SDR's raw sampling freq
#RAW_FS = 200e3					# SDR's raw sampling freq
RAW_FS = 200000                # the queue size is selected so that no more than 1 packet may reside within a queue item
#RAW_FS = 500000           # value should be kept <= expected length, so that we don't miss empty space

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
PREAMBLE_SZ = 1*N*UPSAMP
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
FS = 200000
#DB_THRESH = -7
DB_THRESH = -11
##########################################################################

def spawn_a_worker(my_channel, input_queue, output_queue):
    ###########################################################################################################
    # worker = YourDemodulatorImplementation(my_channel, input_queue, output_queue)
    # worker.start_consuming()
    ###########################################################################################################
    #css_demodulator = CssDemod(N, UPSAMP,PREAMBLE_SZ,END_DELIMITER,DB_THRESH);
    #css_demodulator = CssDemod(N, UPSAMP,PREAMBLE_SZ,END_DELIMITER,DB_THRESH, symbols,PAYLOAD_LEN,NUMPKTS,SF);
    css_demodulator = CssDemod(N, UPSAMP,PREAMBLE_SZ,END_DELIMITER,DB_THRESH, symbols,PAYLOAD_LEN,NUMPKTS,SF,BW,FS);
    outfile = ('../simpleTX_sim/'+exp_root_folder + '/' + exp_folder + '/' + 'error_out')
    css_demodulator.setErrorMeasurementFile(outfile)
    print("Started")
    max_queue_cnt = 10
    while (True):
        if (input_queue.qsize() > 0): 
            queue = input_queue.get()[0]            
            output = []
            css_demodulator.css_demod(my_channel, queue, output)     
            if (len(output) >= 1):
                #print(output)
                print("====")
    print("Done.")
        
def IQ_SOURCE(chan, chan_num):
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
                    chan.put((final, time.time()))
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
        multiprocessing.Process(target=IQ_SOURCE, args=(channel_streams[i], LORA_CHANNELS[i])).start()

    time.sleep(7260)