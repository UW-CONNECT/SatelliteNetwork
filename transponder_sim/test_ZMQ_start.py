'''
If the program is ended pre-maturely, just restart GNURADIO. ZMQ start-of-sync is not the best.
'''
import socket
import numpy as np
import subprocess
import shlex
import zmq

import time
import threading 
# FILE_PATH = r'J:\schellberg\indoor_exp_feb_2024\doppler_simulation_files\GNURADIO_linear_8k_sep'
FILE_PATH = r'..\..\doppler_simulation_files\GNURADIO_linear_8k_sep'
context = zmq.Context() 
socket = context.socket(zmq.PUB) 
socket.bind("tcp://127.0.0.1:4444")


with open(FILE_PATH, "rb") as f:
    nbytes= int((32 /8) *2 * 200000)
    print(nbytes)
    while (byte := f.read(nbytes)):
        # Do stuff with byte.
    
        socket.send(byte)
        byte = np.frombuffer(byte, dtype=np.float32)
        print("Sent byte")
        time.sleep(.1)
print("Sent doppler info!")