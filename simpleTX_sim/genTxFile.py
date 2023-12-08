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
preamble = [1, 1,1, 1,1, 1,1, 1]
#end_delimeter = [1,0,5,0,0,10,1,1] 
end_delimeter = [1] 
css_modulator = CssMod(N, UPSAMP, preamble, end_delimeter) 
output = css_modulator.symbol2packet(symbols)

bin_dat = np.float32(css_modulator.ang2bin(output))
bin_dat = bin_dat.tobytes()

print(np.concatenate((preamble, [49],np.squeeze(symbols), end_delimeter)))

file = open('binaryTXdata_2', 'bw')
file.write(bin_dat)
file.close()
