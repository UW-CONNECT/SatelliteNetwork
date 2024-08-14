'''
If the program is ended pre-maturely, just restart GNURADIO. ZMQ start-of-sync is not the best.
'''
import socket
import numpy as np
# import subprocess
# import shlex
import zmq
import time
import threading 
import math
import system 

os.system('python3 doppler_transponder_ZMQ_noGUI.py &')
FS = 200e3 
doppler_slope=50

# FILE_PATH = r'J:\schellberg\indoor_exp_feb_2024\doppler_simulation_files\GNURADIO_linear_8k_sep'
# FILE_PATH = r'..\..\doppler_simulation_files\GNURADIO_linear_8k_sep'
context = zmq.Context() 
socket = context.socket(zmq.PUB) 
socket.bind("tcp://127.0.0.1:4444")

start_time = time.time()
sweep_time =  2000 # time in seconds 
# with open(FILE_PATH, "rb") as f:
# nbytes= int((32 /8) *2 * 200000)
nsamples = 2*1e6 
# print(nbytes)

t_end = 0 
while (t_end < sweep_time):
    time_vec = t_end + np.arange(0, nsamples) / FS
    t_end = time_vec[-1] + 1/FS
    output = np.exp(1j*(2*math.pi * doppler_slope * (time_vec**2) / 2))
    byte = np.zeros((2*len(output)))
    byte[0::2] = np.real(output)
    byte[1::2] = np.imag(output)
    byte = np.float32(byte)
    #Do stuff with byte.
    socket.send(byte.tobytes())
    # byte = np.frombuffer(byte, dtype=np.float32)
    print("Sent byte")
    # time.sleep(.1)
# print("Sent doppler info!")
print("Finished transmitting doppler for the sweep time.")