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
    def __init__(self, N, SF, BW, FS, preamble, end_delimeter):
        self.N = N 
        
        self.SF = SF 
        
        self.BW = BW 
        
        self.FS = FS 
        
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
       
        #pkt = self.sym_to_data_ang(np.concatenate((self.PREAMBLE , [pkt_length_2], [pkt_length_1], symbols , self.END_DELIMETER)))
        pkt = self.sym_2_css(np.concatenate((self.PREAMBLE , [pkt_length_2], [pkt_length_1], symbols , self.END_DELIMETER)),  self.N, self.SF, self.BW, self.FS)
        print("pkt shape",pkt.shape)
        return pkt   
    
    

    def sym_2_css(self, symbol, N, SF, BW, Fs): 
        '''
        sym_2_css :: Modules a symbol between (0,N-1) to a chirp  
        
        N = 2^SF 
        SF 
        BW = Bandwidth sweep of the chirp -BW/2 to BW/2 
        Fs = Sampling rate
        '''
        sym_out = []
        
        # Fs/Bw is upsamp...
        spsym = int(Fs/BW*N) # symbols defined by their freq offset at the start 
        
        T = spsym/Fs
        k = BW/(T) 
        f_start = -BW/2 
        pi = math.pi 
    
        for sym in symbol: 
            st_offset = int(sym)
            t_arr = np.arange(0, spsym*(N-st_offset)/N)/Fs
            cs = np.exp(1j * 2 * pi *(t_arr*(f_start + k*T*st_offset/N + 0.5*k*t_arr)))   # don't forget residuals
            
            if (len(t_arr) >0):
                ang_start = np.angle(cs[len(t_arr)-1]) # we need to keep track of the phase 
            else: 
                ang_start = 0
                    
            #t_arr = np.arange(0, spsym*st_offset/N - 1)/Fs # [TODO] fix this ... 
            t_arr = np.arange(0, spsym*st_offset/N )/Fs 
            ce = np.exp(1j*(ang_start + 2*pi*(t_arr*(f_start +0.5*k*t_arr)) ))
            out_sig = np.concatenate((cs, ce))
            #sym_out.append(out_sig)
            sym_out=np.append(sym_out,out_sig)
            '''
            plt.figure(1)
            plt.specgram(out_sig) 
            plt.title('out')
            plt.show()
            '''
        return sym_out
    
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