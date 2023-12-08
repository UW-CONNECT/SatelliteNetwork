'''
css_mod.py :: Jackie Schellberg :: November 1, 2023 

Code related to modulating CSS signals 

Current implementation assumes that SF/sampling rate is agreed upon before hand between Tx and Rx code.

'''

from matplotlib import pyplot as plt
import time
import math
import numpy as np
from scipy import signal, fft

import scipy.io
class CssMod:
    def __init__(self, N, UPSAMP, preamble, end_delimeter):
        self.N = N 
        
        self.UPSAMP = UPSAMP 
        
        self.PREAMBLE =  self.sym_to_data_ang(preamble)
        
        self.END_DELIMETER =  self.sym_to_data_ang(end_delimeter)
      
    def symbol2packet(self, symbols): 
        '''
        Converts symbols to binary and appends preamble, packet length, data, and ending delimeter. 
        '''
        output_packet_ang = [] 
        
        pkt_length = len(symbols)       
        if pkt_length > self.N: 
            print("Payload is too big, max size is: ", self.N)
            return 
        
        syms = np.concatenate((self.PREAMBLE , self.sym_to_data_ang([pkt_length]), \
            self.sym_to_data_ang(symbols) , self.END_DELIMETER))
        
        
        return syms   
        #return output_packet_ang
        
    def sym_to_data_ang(self, symbol):
            '''
            Via https://github.com/mananmishra11/open-lora/blob/main/std_lora/sym_to_data_ang.py
            '''
            data = []
            out_data = []
            #accumulator = 0
            pi = math.pi

            for j in symbol:
                accumulator = 0
                phase = -pi + (j-1)*(2*pi/(self.N))
                temp = np.zeros((self.N, 1), complex)
                for i in range(0, self.N):
                    accumulator = accumulator + phase
                    polar_radius = 1
                    
                    x = polar_radius * math.cos(accumulator)
                    y = polar_radius * math.sin(accumulator)
                    temp[i] = complex(x, y)
                    phase = phase + (2*pi/self.N)
                    
                data = temp
                
                data = np.squeeze(data)
                
                # downsample the signal by adding zeros in the fft 
                if(self.UPSAMP != 0):
                    data_fft = fft.fft(data)
                    nz_h_1 = data_fft[:int(len(data_fft)/2)] 
                    zeroes_h = np.squeeze(np.zeros(((self.UPSAMP-1)*len(data_fft), 1), complex))
                    nz_h_2 = data_fft[int(len(data_fft)/2):]
                    parts = np.concatenate((nz_h_1, zeroes_h , nz_h_2))
                    data_upsamp = fft.ifft(parts)           
                    data = data_upsamp
                     
                out_data = np.append(out_data, data)
            return out_data
            
    def ang2bin(self, data):
        '''
        Convert complex data to binary for GNUradio 
        '''
        Rx_buff = np.zeros((2*len(data)))
        Rx_buff[0::2] = np.real(data)
        Rx_buff[1::2] = np.imag(data)
        Rx_buff = np.append(Rx_buff, np.zeros((len(Rx_buff))))
        
        return Rx_buff