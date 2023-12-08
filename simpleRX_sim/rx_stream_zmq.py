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
LOGGER_LOC = 'DATA_LOG.txt'
#RAW_FS = 450e3					# SDR's raw sampling freq
#RAW_FS = 200e3					# SDR's raw sampling freq
#RAW_FS = 200000                 # the queue size is selected so that no more than 1 packet may reside within a queue item
RAW_FS = 200000 
LORA_BW = 125e3        # LoRa bandwidth
LORA_CHANNELS = [1]  # channels to process

UPSAMPLE_FACTOR = 4             		# factor to downsample to
#OVERLAP = 10 * UPSAMPLE_FACTOR#int(10 * RAW_FS / BW)
OVERLAP = 0

# CSS Specifics to be known ahead of time
SF = 7 
N = 2**SF
UPSAMP = 10;

##########################################################################

def spawn_a_worker(my_channel, input_queue, output_queue):
    ###########################################################################################################
    # worker = YourDemodulatorImplementation(my_channel, input_queue, output_queue)
    # worker.start_consuming()
    ###########################################################################################################
    css_demodulator = CssDemod(N, UPSAMP);
    
    print("Started")
    max_queue_cnt = 10
    #while (output_queue.qsize() < max_queue_cnt):
    while (True):
        if (input_queue.qsize() > 0): 
            queue = input_queue.get()[0]            
            output = []
            css_demodulator.css_demod(my_channel, queue, output)         
            #output_queue.put(output)
            '''
            plt.figure(1)
            xpts = range(0,len(queue))
            plt.plot(xpts, 20*np.log10(queue))
            plt.show()
              '''  
            if (len(output) > 1):
                print(output)
                '''
                plt.figure(1)
                xpts = range(0,len(queue))
                plt.plot(xpts, 10*np.log10(queue))
                
                plt.figure(2)
                xpts = range(0,len(output))
                plt.plot(xpts, output)
                
                plt.show()
                '''
    '''
    for i in range (0, max_queue_cnt):
        print("Displaying the queue items")
        output2 = output_queue.get()
        print(len(output2))
        plt.figure(1)
        xpts = range(0,len(output2))
        plt.plot(xpts, output2)
        ax = plt.gca()
        ax.set_ylim([0, 128])
        plt.title("Demodulated symbols") 
        plt.show()      
    '''
    print("Done.")
        
def IQ_SOURCE(chan, chan_num):
    #from client_config import RAW_FS, OVERLAP
    t2 = '127.0.0.1'
    if chan_num == 1:
        port = 4900
    '''
    server_socket = socket.socket(socket.AF_INET, socket.SOCK_DGRAM)
    server_socket.bind((t2, port))
    server_socket.setsockopt(socket.SOL_SOCKET, socket.SO_RCVBUF, 2 ** 20)
    '''
    context = zmq.Context() 
    socket = context.socket(zmq.SUB)
    socket.connect("tcp://127.0.0.1:55555")
    socket.setsockopt(zmq.SUBSCRIBE, b'')
    
    recver = list()
    counter = 0
    tmp = np.array([])
    atime = time.time()

    while True:
        #message, address = server_socket.recvfrom(2 ** 15)
        if socket.poll(10) != 0:
            message = socket.recv() 
            '''
            data = np.frombuffer(message, dtype=np.complex64, count=-1) 
            plt.figure(1)
            xpts = range(0,len(data))
            plt.plot(xpts, abs(data))
            ## [JMS] : Come back here to make sure we are decoding the data properly
            plt.show()
            '''
            
            #print("Length of the socket msg: ", len(message))
            if len(message) > 0:
                #print(len(message))
                counter += len(message)
                data = np.frombuffer(message, dtype=np.complex64, count=-1) 
                #recver.append(message)
                recver.append(data)
                
                
                if counter > int(RAW_FS + OVERLAP * 2) * 8:                 
                    #print("recvelen",len(recver))
                    plt.show()
                    final = np.concatenate(recver).ravel()
                    #final = np.array(recver).view(dtype=np.complex64)       # flatten
                    #final = np.array(recver).view(dtype=np.float32)       # flatten
                    #final = final[::2] + 1j*final[1::2] 
                    #final = np.frombuffer(np.array(recver), dtype=np.complex64, count=-1) 
                    final = np.concatenate((tmp, final))                    # combine with previous
                    xl =len(tmp)
                    tmp = final[int(RAW_FS):]                               # get overlap
                    #print("TMP LENGTH...",len(tmp))
                    #final = final[0:int(RAW_FS + OVERLAP * 2)]              # final chunk
                    final = final[0:int(RAW_FS)]
                    chan.put((final, time.time()))
                    #print(f"Received a chunk at time {time.time()}\n")      # Comment this print later
                    recver = list()
                    counter = len(tmp) * 8
                    atime += 1
                    
                    '''
                    plt.figure(1)
                    xpts = range(0,len(final))
                    plt.plot(xpts, abs(final))
                    plt.axhline(y =0, color = 'b')
                    plt.axvline(x=xl, color = 'b')
                    
                    cx = range(xl,len(final),4096)
                    for cc in cx:
                        plt.axvline(x=cc, color = 'b',ls=':')
                    plt.show()
                    '''
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