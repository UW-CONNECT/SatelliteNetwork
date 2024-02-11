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
        
        self.PREAMBLE =  preamble
        
        self.END_DELIMETER =  end_delimeter
      
    def symbol2packet(self, symbols): 
        '''
        Converts symbols to binary and appends preamble, packet length, data, and ending delimeter. 
        '''
        output_packet_ang = [] 
        
        pkt_length_2 = len(symbols) % self.N 
        pkt_length_1 = math.floor(len(symbols) / self.N)
        
        if pkt_length_2 > self.N: 
            print("Payload is too big, max size is: ", self.N)
            #return 
       
        pkt = self.sym_to_data_ang(np.concatenate((self.PREAMBLE , [pkt_length_2], [pkt_length_1], symbols , self.END_DELIMETER)))
        
        return pkt   
        
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
                #samps = self.N * self.UPSAMP 
                
                samps = self.N 
                phase = -pi + (j-1)*(2*pi/samps)
                temp = np.zeros((samps, 1), complex)
                #for i in range(0, self.N):
                
                for i in range(0, samps):
                    accumulator = accumulator + phase
                    polar_radius = 1
                    
                    x = polar_radius * math.cos(accumulator)
                    y = polar_radius * math.sin(accumulator)
                    temp[i] = complex(x, y)
                    #phase = phase + (2*pi/self.N)
                    phase = phase + (2*pi/samps)
                    
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
                    #data = data_upsamp * self.UPSAMP
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
        #Rx_buff = np.append(np.zeros((len(Rx_buff))),Rx_buff)
        
        return Rx_buff
        
    def ang2bin_nopad(self, data):
        '''
        Convert complex data to binary for GNUradio 
        '''
        Rx_buff = np.zeros((2*len(data)))
        Rx_buff[0::2] = np.real(data)
        Rx_buff[1::2] = np.imag(data)
       
        #Rx_buff = np.append(np.zeros((len(Rx_buff))),Rx_buff)
        
        return Rx_buff