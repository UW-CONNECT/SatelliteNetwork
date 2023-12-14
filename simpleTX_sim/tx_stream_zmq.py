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

SF = 7 
N = 2**SF
UPSAMP = 10;
symbols = np.concatenate((np.ones(8), range(1,50)))

preamble = [1,1] 
end_delimeter = [3,3,3,3]
css_modulator = CssMod(N, UPSAMP,preamble, end_delimeter) 
output = css_modulator.symbol2packet(symbols)
print(len(output))
bin_dat = np.float32(css_modulator.ang2bin(output))
bin_dat = bin_dat.tobytes()

# For testing, continuously transmit these symbols 
t = time.time()
'''
#socket.send(bin_dat)
file = open('outdat_2', 'bw')

bufsize=10000
endloc = bufsize 
num_pkts = math.ceil(len(bin_dat)/bufsize)
for i in range(0, num_pkts):
    if (i == (num_pkts-1)): 
        endloc = len(bin_dat[i*bufsize:])
        
    #socket.send(bin_dat[i*bufsize:i*bufsize+endloc])
    file.write(bin_dat[i*bufsize:i*bufsize+endloc])
    time.sleep(.01)
  
print(len(bin_dat))

file.close()
'''
time.sleep(1)
#while(True):
socket.send(bin_dat)
#    time.sleep(72960/200000)


print(time.time() - t)
print("Packets sent.")
#file.close()
#time.sleep(3)

# test to make sure modulation works by demodulating w/o any channel effects  
'''
plt.figure(1)
xpts = range(0,len(output))
plt.plot(xpts, np.real(output))
plt.show()  

import sys
sys.path.append('../simpleRX_sim')
from css_demod import CssDemod
SF = 7 
N = 2**SF
UPSAMP = 10;
css_demodulator = CssDemod(N, UPSAMP);

out_sim = []
css_demodulator.css_demod([], output, out_sim) 

print(len(out_sim))

out_sim = np.squeeze(out_sim)
plt.figure(1)
xpts = range(0,len(out_sim))
plt.plot(xpts, out_sim)
plt.show()  '''

# we can send out samples ASAP and let GNUradio throttle to our desired sampling rate (200k)



