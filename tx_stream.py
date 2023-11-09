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

# offsets that need to be considered with payload size
IP_HEADER_SIZE = 1 
IP_PACKET_HEADER = 20
UDP_DGRAM_HEADER_SIZE = 8 

UDP_IP = "127.0.0.1"
UDP_PORT = 4900
sock = socket.socket(socket.AF_INET,
                     socket.SOCK_DGRAM) 
                     
#bufsize = sock.getsockopt(socket.SOL_SOCKET, socket.SO_SNDBUF)
#sock.setsockopt(socket.SOL_SOCKET, socket.SO_SNDBUF, 1500)   
sock.setsockopt(socket.SOL_SOCKET, socket.SO_SNDBUF, 528)      
bufsize = sock.getsockopt(socket.SOL_SOCKET, socket.SO_SNDBUF)
#bufsize = bufsize-8

#time.sleep(10)
bufsize = bufsize - (IP_HEADER_SIZE+IP_PACKET_HEADER+UDP_DGRAM_HEADER_SIZE) 

print(bufsize)
SF = 7 
N = 2**SF
UPSAMP = 10;
#symbols = [0, 0, 0, 0, 0, 0, 0, 0, 1, 2, 3 , 4, 5, 6, 7, 8 ]
symbols = np.ones(50)*2

css_modulator = CssMod(N, UPSAMP) 
output = css_modulator.symbol2packet(symbols)
bin_dat = np.float32(css_modulator.ang2bin(output))
bin_dat = bin_dat.tobytes()

# For testing, continuously transmit these symbols 
#while(True):
#for i in range(0, math.ceil(len(bin_dat)/bufsize)):
t = time.time()
endloc = bufsize 
num_pkts = math.ceil(len(bin_dat)/bufsize)
print(len(bin_dat))

# use this for comparison
#file = open('outdat_2', 'bw')
#file.write(bin_dat)
#file.close()

for i in range(0, num_pkts):
    #print(len(bin_dat[i*bufsize:i*bufsize+bufsize]))
    
    # what's a good way to verify the speed?
    
    # convert the floats to bytes       
    #sock.sendto(bin_dat[i*bufsize:i*bufsize+bufsize], (UDP_IP, UDP_PORT))
    '''
    print(i*bufsize)
    print(i*bufsize+bufsize)
    print('End')
    '''
    #print(i)
    if (i == (num_pkts-1)): 
        endloc = len(bin_dat[i*bufsize:])
        
    sock.sendto(bin_dat[i*bufsize:i*bufsize+endloc], (UDP_IP, UDP_PORT))
    #file.write(bin_dat[i*bufsize:i*bufsize+endloc])
    time.sleep(.01)

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



