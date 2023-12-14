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

''' load ref packet for BER testing '''
gndtruth = [16, 79, 98, 99, 49, 113, 94, 91, 109, 55, 18, 89, 49, 49, 49, 108, 84, 84, 15, 12, 81, 77, 103, 50, 19, 80, 13, 73, 82, 6, 125, 121, 119, 59, 13, 37, 86, 87, 69, 87, 71, 21, 44, 62, 114, 26, 52, 15, 2, 117, 26, 43, 98, 61, 18, 32, 119, 46, 116, 122, 119, 8, 79, 37, 104, 12, 47, 28, 43, 80, 60, 14, 74, 106, 68, 38, 69, 92, 89, 10, 66, 27, 69, 11, 42, 8, 83, 122, 67, 126, 125, 38, 27, 113, 66, 102, 61, 124, 109, 84, 93, 5, 33, 37, 124, 100, 8, 76, 105, 46, 45, 118, 23, 96, 107, 12, 95, 85, 59, 89, 118, 116, 120, 60, 32, 64, 122, 33, 19, 95, 62, 63, 13, 102, 1, 102, 120, 99, 42, 44, 74, 99, 39, 117, 88, 26, 2, 67, 14, 37, 15, 123, 42, 10, 21, 125, 45, 23, 90, 74, 112, 77, 89, 107, 125, 94, 68, 73, 86, 124, 121, 24, 107, 95, 12, 83, 75, 60, 73, 124, 56, 2, 10, 75, 83, 53, 72, 120, 107, 38, 8, 77, 40, 54, 116, 50, 120, 65, 38, 120, 74, 27, 81, 3, 108, 59, 33, 95, 77, 33, 97, 41, 109, 16, 102, 15, 42, 84, 52, 91, 118, 34, 48, 33, 109, 20, 92, 23, 10, 94, 13, 126, 15, 84, 7, 114, 34, 99, 78, 39, 103, 101, 41, 18, 120, 39, 34, 46, 13, 17, 37, 72, 89, 8, 124, 10, 74, 108, 29, 125, 26, 85, 74, 54, 111, 96, 20, 43, 85, 51, 101, 83, 40, 22, 66, 51, 11, 83, 67, 23, 111, 35, 68, 42, 48, 16, 28, 49, 20, 5, 15, 13, 40, 27, 7, 32, 22, 55, 41, 126, 82, 28, 78, 19, 105, 84, 121, 28, 47, 39, 91, 37, 31, 4, 99, 47, 1, 52, 83, 63, 4, 37, 43, 2, 109, 6, 72, 18, 110, 110, 11, 85, 4, 5, 76, 10, 117, 85, 104, 29, 58, 78, 124, 27, 113, 9, 28, 5, 2, 19, 38, 25, 16, 101, 64, 48, 111, 75, 69, 109, 57, 114, 110, 25, 58, 109, 40, 27, 9, 48, 41, 70, 53, 67, 3, 124, 93, 41, 29, 90, 86, 79, 117, 73, 120, 63, 29, 82, 49, 28, 121, 62, 28, 127, 20, 67, 125, 121, 6, 104, 1, 101, 15, 105, 101, 114, 89, 13, 18, 43, 96, 116, 106, 99, 79, 69, 125, 102, 58, 100, 53, 39, 81, 50, 57, 51, 106, 67, 48, 106, 95, 40, 29, 110, 23, 91, 70, 5, 29, 92, 65, 47, 100, 14, 111, 20, 37, 122, 24, 84, 23, 44, 72, 14, 44, 27, 78, 31, 108, 79, 87, 13, 40, 63, 24, 108, 100, 49, 28, 100, 73, 42, 72, 110, 2, 65, 75, 120, 68, 65, 30, 121, 1, 54, 23, 45, 105, 111, 5, 21, 32, 87, 15, 50, 29, 62, 107, 73, 118, 85, 103, 37, 44, 111, 52, 115, 79, 33, 28, 95, 40, 92, 95, 45, 57, 12, 89, 46, 69, 121, 76, 94, 82, 10, 46, 95, 127, 49, 71, 103, 35, 6, 10, 108, 58, 53, 66, 98, 15, 42, 88, 83, 7, 95, 66, 62, 26, 13, 53, 27, 119, 123, 31, 124, 99, 119, 3, 16, 105, 119, 15, 60, 68, 78, 127, 97, 7, 110, 57, 38, 19, 118, 105, 102, 34, 16, 92, 24, 8, 58, 42, 83, 5, 63, 32, 110, 15, 84, 44, 35, 2, 8, 40, 120, 18, 68, 71, 73, 34, 124, 98, 91, 70, 49, 22, 11, 97, 99, 78, 81, 98, 55, 33, 112, 115, 55, 8, 88, 53, 91, 88, 14, 81, 108, 11, 63, 52, 11, 48, 23, 50, 85, 9, 73, 13, 14, 114, 122, 55, 42, 60, 98, 93, 15, 99, 114, 119, 2, 54, 84, 13, 14, 28, 21, 4, 3, 113, 22, 37, 81, 95, 112, 116, 44, 7, 33, 8, 73, 71, 72, 40, 78, 74, 37, 43, 76, 101, 42, 76, 84, 25, 27, 45, 21, 94, 49, 68, 33, 32, 125, 32, 63, 46, 126, 71, 124, 80, 15, 63, 31, 106, 87, 23, 95, 87, 27, 9, 85, 68, 110, 28, 30, 108, 60, 56, 87, 8, 71, 110, 34, 16, 66, 36, 51, 78, 22, 93, 118, 51, 98, 82, 67, 51, 2, 105, 30, 7, 119, 56, 19, 68, 92, 6, 106, 76, 11, 15, 25, 103, 14, 9, 95, 23, 105, 18, 11, 97, 127, 53, 41, 36, 31, 98, 121, 59, 55, 82, 45, 81, 34, 16, 50, 112, 8, 98, 41, 52, 22, 94, 57, 109, 108, 84, 40, 21, 47, 51, 39, 92, 39, 26, 13, 79, 14, 9, 96, 9, 7, 58, 28, 91, 79, 54, 45, 113, 8, 24, 124, 4, 5, 38, 6, 38, 52, 87, 102, 40, 59, 81, 97, 88, 65, 61, 76, 105, 82, 38, 123, 109, 44, 28, 120, 50, 57, 125, 124, 79, 107, 103, 17, 100, 113, 26, 82, 57, 112, 9, 51, 38, 102, 56, 5, 19, 59, 74, 111, 9, 18, 52, 37, 37, 51, 102, 7, 98, 118, 115, 53, 42, 67, 84, 48, 114, 26, 73, 40, 20, 9, 33, 53, 67, 48, 76, 8, 6, 78, 36, 24, 123, 1, 80, 58, 28, 111, 47, 79, 80, 91, 35, 104, 106, 55, 55, 119, 77, 109, 33, 110, 47, 73, 113, 94, 78, 36, 37, 58, 46, 51, 48, 25, 68, 75, 127, 37, 103, 58, 75, 22, 21, 28, 66, 39, 4, 43, 46, 105, 105, 98, 103, 8, 72, 33, 42, 4, 71, 54, 77, 45, 73, 13, 16, 24, 74, 104, 55, 86, 90, 11, 31, 5, 26, 112, 107, 34, 100, 27, 2, 9, 12, 2, 52, 59, 69, 36, 105, 57, 108, 5, 69, 114, 126, 99, 26, 8, 44, 110, 117, 34, 115, 106, 37, 1, 126, 67, 26, 123, 22, 36, 49, 40]

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
                print(np.array(output) - np.array(gndtruth))
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