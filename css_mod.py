'''
css_mod.py :: Jackie Schellberg :: November 1, 2023 

Code related to modulating CSS signals 

Current implementation assumes that SF/sampling rate is agreed upon before hand.

'''

from matplotlib import pyplot as plt
import time
import math
import numpy as np
from scipy import signal, fft

import scipy.io
class CssMod:
    def __init__(self, N, UPSAMP):
        self.N = N 
        
        self.UPSAMP = UPSAMP 
        
        self.PREAMBLE =  self.sym_to_data_ang([0,0,0,0,0,0,0,0])
        
        self.END_DELIMETER =  self.sym_to_data_ang([0,0,0,0,0,0,0,0])
      
    def symbol2packet(self, symbols): 
        '''
        Converts symbols to binary and appends preamble, packet length, data, and ending delimeter. 
        '''
        output_packet_ang = [] 
        #pkt_length = len(symbol) * self.N * self.UPSAMP
        pkt_length = len(symbols)       
        
        
        syms = self.sym_to_data_ang(symbols)
        
        '''
        # need to check if packet length is within the modulation order, 
        # or we will need multiple chirps to encode it 
        if (pkt_length < self.N): 
            pkt_size = self.sym_to_data_ang([pkt_length])
            output_packet_ang = np.concatenate([self.PREAMBLE,self.PREAMBLE,pkt_size,syms,self.END_DELIMETER])
        else: 
            print("MO not supported. ")
            output_packet_ang = np.concatenate([self.PREAMBLE,self.PREAMBLE,syms,self.END_DELIMETER])
            #return 0 
        '''    
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
            #out_data.flatten()
            return out_data
            
    def ang2bin(self, data):
        '''
        Convert complex data to binary for GNUradio 
        '''
        #Rx_buff = np.zeros((2*len(data)),complex)
        Rx_buff = np.zeros((2*len(data)))
        Rx_buff[0::2] = np.real(data)
        Rx_buff[1::2] = np.imag(data)
        #Rx_buff = np.append(Rx_buff, np.zeros((len(Rx_buff)),complex))
        Rx_buff = np.append(Rx_buff, np.zeros((len(Rx_buff))))
        
        return Rx_buff