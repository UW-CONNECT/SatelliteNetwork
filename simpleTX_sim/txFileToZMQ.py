
import multiprocessing
from multiprocessing import Manager

from threading import Thread
from queue import Queue
import zmq
import socket
import threading 
import time
import numpy as np
import subprocess
EXP_ROOT_FOLDER ='../../experiment_data_synced'
RAW_DATA_FILE_NAME = 'trial1'
TRIAL_FOLDER = '0'

f = open('../trial_under_test.txt')
TRIAL_FOLDER = f.read()
f.close()

 # = 
exp_folder = EXP_ROOT_FOLDER + '\\' + TRIAL_FOLDER + '\\' + RAW_DATA_FILE_NAME
def send_tx_zmq_data(fp): 
    '''
    Send doppler data over ZMQ to the doppler transponder. 
    '''
    FILE_PATH = fp
    context = zmq.Context() 
    zsocket = context.socket(zmq.PUB) 
    zsocket.bind("tcp://127.0.0.1:4444")
    experiment_start_time = time.time()
    
    experiment_time_path= EXP_ROOT_FOLDER + '\\' + str(TRIAL_FOLDER) + '\\' + "posix_start_time"
    f = open(experiment_time_path, "w")
    f.write(str(experiment_start_time))
    f.close()
    
    print("Started sending TX data over TCP. Start time: ", experiment_start_time)
    with open(FILE_PATH, "rb") as f:
        nbytes= 200000
        print(nbytes)
        while (byte := f.read(nbytes)):
            # Do stuff with byte.        
            zsocket.send(byte)
            byte = np.frombuffer(byte, dtype=np.float32)
            #print("Sent byte")
            time.sleep(.1)
    print("Done sending the data over ZMQ")
    
if __name__ == "__main__":
    
    
    # f = open(experiment_time_path, "w")
    # f.write(str(experiment_start_time))
    # f.close()

    send_tx_zmq_data(exp_folder)
    