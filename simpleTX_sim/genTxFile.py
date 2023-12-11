'''
tx_stream :: Jackie Schellberg :: Last updated Nov 1, 2023 

For streaming a single set of 'symbols' over UDP to GNUradio 
'''

import socket 
from css_mod import CssMod
from matplotlib import pyplot as plt
import numpy as np
import time
import math
import zmq

context = zmq.Context() 
socket = context.socket(zmq.PUB) 
socket.bind("tcp://127.0.0.1:4444")

#SF = 7 
SF = 7 
N = 2**SF
UPSAMP = 10;
#symbols = np.concatenate((np.ones(8), range(1,50)))
symbols =  range(1,50)

#preamble = [0,0,0,0,0,0,0,0]
#end_delimeter = [0,0,0,0,0,0,0,0] 
preamble = [1, 1]
#end_delimeter = [1,0,5,0,0,10,1,1] 
end_delimeter = [3,3,3,3] 
css_modulator = CssMod(N, UPSAMP, preamble, end_delimeter) 
output = css_modulator.symbol2packet(symbols)

bin_dat = np.float32(css_modulator.ang2bin(output))
#print(bin_dat)


bsl = len(bin_dat)*2
bin_dat = np.concatenate((bin_dat,bin_dat,bin_dat,bin_dat,bin_dat))
#print(bin_dat.shape)
#bin_dat = np.concatenate((np.zeros((bsl)), bin_dat))
#print(bin_dat.shape)
#print(len(bin_dat))
bin_dat = np.append(np.float32(np.zeros((bsl))),bin_dat)
bin_dat = bin_dat.tobytes()

#print(np.concatenate((preamble, [49],np.squeeze(symbols), end_delimeter)))



file = open('5longpkts_2oneP_3333', 'bw')
file.write(bin_dat)
file.close()
