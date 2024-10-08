'''
tx_stream :: Jackie Schellberg :: Last updated Nov 1, 2023 

For streaming a single set of 'symbols' over UDP to GNUradio 
'''
from css_mod import CssMod
from matplotlib import pyplot as plt
import numpy as np
import time
import math
import zmq
import random

#SF = 7 
SF = 9
N = 2**SF
#UPSAMP = 10;
UPSAMP = 10;
#symbols = np.concatenate((np.ones(8), range(1,50)))
#symbols =  range(1,50)
#symbols = random.sample(range(1, N), 1000)

fileout_name = "test_hamming_encode"

symbols = []
#for ss in range(0,1000):
#    symbols.append(random.randint(1,N-1))
# symbols = np.ones(100)*5
symbols = [5]

#preamble = [0,0,0,0,0,0,0,0]
#end_delimeter = [0,0,0,0,0,0,0,0] 
#preamble = [1, 0,0,1,1,0,0,1]
preamble = [1]
#preamble = [1,1,1,1]
#end_delimeter = [1,0,5,0,0,10,1,1] 
end_delimeter = [3,3,3,3] 
CR = 3  # CR=0 means no hamming codes 
BW = 20000
FS = 200000
css_modulator = CssMod(N, SF, BW, FS, preamble, end_delimeter, CR) 
output = css_modulator.symbol2packet(symbols)

bin_dat = np.float32(css_modulator.ang2bin(output))
#print(bin_dat)


bsl = len(bin_dat)
#bin_dat = np.concatenate((bin_dat,bin_dat,bin_dat,bin_dat,bin_dat))

#print(bin_dat.shape)
#bin_dat = np.concatenate((np.zeros((bsl)), bin_dat))
#print(bin_dat.shape)
#print(len(bin_dat))
buffer_dat = N*UPSAMP*50*32
bin_dat = np.append(np.float32(np.zeros((buffer_dat))),bin_dat)
#bin_dat = np.append(np.float32(np.zeros((bsl))),bin_dat)
bin_dat = bin_dat.tobytes()

print(len(bin_dat))
#print(np.concatenate((preamble, [49],np.squeeze(symbols), end_delimeter)))

file = open(fileout_name, 'bw')
file.write(bin_dat)
file.close()

file = open(fileout_name + '_ref','w')
file.write(f"{symbols}")
file.close()
