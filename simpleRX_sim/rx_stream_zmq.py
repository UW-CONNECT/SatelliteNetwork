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
from css_demod import CssDemod

import zmq

from matplotlib import pyplot as plt

############################### Variables ###############################
#RAW_FS = 450e3					# SDR's raw sampling freq
#RAW_FS = 200e3					# SDR's raw sampling freq
#RAW_FS = 200000                # the queue size is selected so that no more than 1 packet may reside within a queue item
RAW_FS = 200000           # value should be kept <= expected length, so that we don't miss empty space

LORA_CHANNELS = [1]  # channels to process

UPSAMPLE_FACTOR = 4             		# factor to downsample to
#OVERLAP = 10 * UPSAMPLE_FACTOR#int(10 * RAW_FS / BW)
OVERLAP = 0

# CSS Specifics to be known ahead of time
SF = 7 
N = 2**SF
UPSAMP = 10;
PREAMBLE_SZ = N*UPSAMP
END_DELIMITER = [3,3,3,3] 

##########################################################################

def spawn_a_worker(my_channel, input_queue, output_queue):
    ###########################################################################################################
    # worker = YourDemodulatorImplementation(my_channel, input_queue, output_queue)
    # worker.start_consuming()
    ###########################################################################################################
    css_demodulator = CssDemod(N, UPSAMP,PREAMBLE_SZ,END_DELIMITER);
    
    print("Started")
    max_queue_cnt = 10
    while (True):
        if (input_queue.qsize() > 0): 
            queue = input_queue.get()[0]            
            output = []
            css_demodulator.css_demod(my_channel, queue, output)     
            if (len(output) >= 1):
                print(output)
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