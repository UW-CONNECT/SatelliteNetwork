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
    def __init__(self, N, SF, BW, FS, preamble, end_delimeter, CR):
        self.N = N 
        
        self.SF = SF 
        
        self.BW = BW 
        
        self.FS = FS 
        
        self.PREAMBLE =  preamble
        
        self.END_DELIMETER =  end_delimeter
        
        self.CR = CR
      
    def symbol2packet(self, symbols): 
        '''
        Converts symbols to binary and appends preamble, packet length, data, and ending delimeter. 
        '''
        output_packet_ang = [] 
        
        symbols_tmp = symbols
        symbols = self.encode(symbols,self.CR)
        
        pkt_length_2 = len(symbols) % self.N 
        pkt_length_1 = math.floor(len(symbols) / self.N)
        
        print(pkt_length_2 + self.N*pkt_length_1, len(symbols))
        # pkt_length_2 = math.ceil(len(symbols) / self.N) 
        # pkt_length_1 = math.floor(len(symbols) / self.N)
        print("num wrong: ",sum(np.subtract(self.decode(symbols,self.CR), symbols_tmp)))
        if pkt_length_2 > self.N: 
            print("Payload is too big, max size is: ", self.N)
            #return 
       
        #pkt = self.sym_to_data_ang(np.concatenate((self.PREAMBLE , [pkt_length_2], [pkt_length_1], symbols , self.END_DELIMETER)))
        
        # payload = np.concatenate((self.PREAMBLE , [pkt_length_2], [pkt_length_1], symbols , self.END_DELIMETER))
        
        sync_word = np.conjugate(self.sym_2_css([0, 0], self.N, self.SF, self.BW, self.FS))
        # payload = np.concatenate((self.PREAMBLE , sync_word, [pkt_length_2], [pkt_length_1], symbols , self.END_DELIMETER))        
       
        # if self.CR >= 1: 
        
        # print("Encoded:",payload)
        # payload[5] = payload[5] + 123
        # print("Error:",payload)
        # payload = self.decode(payload)
        
        # print(payload)
        pkt = np.concatenate(( self.sym_2_css(self.PREAMBLE,  self.N, self.SF, self.BW, self.FS), sync_word, self.sym_2_css(np.concatenate(([pkt_length_2], [pkt_length_1], symbols , self.END_DELIMETER)),  self.N, self.SF, self.BW, self.FS)))
        
        # plt.figure(1)
        # plt.specgram(pkt)
        # plt.show()
        
        
        # print("pkt shape",pkt.shape)
        return pkt 
        
    def decode(self, payload, CR): 
        '''
        Unpack the nibbles from the symbols and decode  
        '''
        decoded_out = []
        bit_offset = 0 
        sf_dx = self.SF # scrolls for 1 to SF to select bits 
        
        # max_bits_per_sym = (self.SF + self.CR) -1
        max_bits_per_sym = (3 + CR) 
        current_bit = max_bits_per_sym 
        outword = []
        codeword = 0
        # print(payload)
        for sym in payload: 
            # print(sym)
            # for bit in range(0,self.SF):
            for bit in range(self.SF-1,-1,-1):
                # codeword.append(int(sym) & 2**current_bit) 
                # codeword = codeword + (int(sym) & 2**current_bit) 
                has_ones= (sym & 2**(bit)) > 0
                codeword = codeword + (has_ones*2**current_bit) 
                
                # if (current_bit >= max_bits_per_sym):
                if (current_bit <= 0):
                    # print(codeword)
                    # outword = self.hamming_decode([codeword], self.CR)
                    # print("Outword:", outword)
                    # decoded_out.append(outword)
                    decoded_out.append(codeword)
                    current_bit = max_bits_per_sym
                    codeword = 0
                else: 
                    current_bit = current_bit - 1
                
        # payload length may not match, check for incomplete 
        if (current_bit != max_bits_per_sym):
            # decoded_out.append(self.hamming_decode([codeword], self.CR))
            decoded_out.append(codeword)          
          
        decoded_out = self.hamming_decode(decoded_out, self.CR)
        # re-pack into SF-bit-wide symbols 
        # now select every 'SF' bits
        nsymbols = int(np.ceil((self.SF + self.CR) * len(payload)) / self.SF)
        symbols = []
        nib_dx = 0
        # nib_offset =  (3 + self.CR) 
        nib_offset =  (3 ) 
        sf_bit = self.SF - 1
        current_symbol = 0
        for nib in decoded_out:
            for bit_ptr in range(nib_offset,-1,-1):
                has_ones= (nib & 2**(bit_ptr)) > 0
                current_symbol = current_symbol + has_ones*2**sf_bit
                
                # current_symbol = 
                if sf_bit <= 0:
                    sf_bit = self.SF - 1
                    symbols.append(current_symbol)
                    current_symbol = 0
                else: 
                    sf_bit = sf_bit - 1
                    
        # payload length may not match, check for incomplete 
        if (sf_bit != self.SF - 1 and current_symbol != 0):
            symbols.append(current_symbol)
        return symbols
        
    def hamming_decode(self, codewords, CR):
        '''
        Hamming decode, where CR is the coding rate. 
        '''
        CR = CR + 4
        nibs = []
        for i in range(0, len(codewords)):
            codeword = codewords[i]
            p1 = self.bit_reduce(codeword,[7,3,2,0])
            p2 = self.bit_reduce(codeword,[6,3,1,0])
            p3 = self.bit_reduce(codeword,[4,2,1,0])
            p4 = self.bit_reduce(codeword,[4,3,2,1,0])
            p5 = self.bit_reduce(codeword,[5,3,2,1])            
            
            # implement this 
            if CR == 5 or CR == 6:
                # detect but don't correct 
                nibs.append(np.mod(codeword, 16)) 
            elif CR == 7 or CR == 8: 
                par = p2*4 + p3*2 + p5 
                # fix the parity and modify par  
                
                if par == 3: 
                    par = 4 
                elif par == 5:
                    par = 8
                elif par == 6:
                    par = 1 
                elif par == 7:
                    par = 2
                else:
                    par = 0
                
                codeword = codeword ^ par 
                nibs.append(np.mod(codeword, 16)) 
                # pass
            else:
                nibs.append(codeword) 
        return nibs 
        
    def hamming_encode(self, nibbles, CR):
        '''
        These are setup to take 8-bit ints. May need to change this later. 
        
        nibbles = 4 bit 
        CR = coding rate (1-4) 
        '''
        n_count = len(nibbles)
        codewords = []
        
        for i in range(0, n_count):
            nib = nibbles[i] 
            p1 = self.bit_reduce(nib, [0, 2, 3])
            p2 = self.bit_reduce(nib, [0, 1, 3])
            p3 = self.bit_reduce(nib, [0, 1, 2])
            p4 = self.bit_reduce(nib, [0, 1, 2, 3])
            p5 = self.bit_reduce(nib, [1, 2, 3])
            
            if (CR == 1):
                codewords.append(p4 * 16 | nib)
            elif (CR == 2):
                codewords.append(p5*32 | p3*16 | nib)
            elif (CR == 3):
                codewords.append(p2*64 | p5*32 | p3*16 | nib)
            elif (CR == 4):
                codewords.append(p1*128 | p2*64 | p5*32 | p3*16 | nib)
            else: 
                # print("Check CR ... [encode]")
                codewords.append(nib)
            
        return np.array(codewords) 
    
    def bit_reduce(self, w,pos):
        '''
        XOR two integers -> Essentially 'and' w with pos, and count the number of 1's mod 2 
        
        0-based indexing for pos[i]
        '''
        b = 0
        #pos = pos - 1
        for i in range(0, len(pos)):
            if ((w & 2**(pos[i])) > 0):
                b = b + 1
        return np.mod(b,2)==1
        
    def encode(self, payload, CR):
        '''
        Prep the payload for hamming encoding, based on the SF. 
        '''
        # convert integers to SF bit numbers, and divide as nibbles 
        payload = np.array(payload,dtype="uint16")
        nmask = int(np.ceil(self.SF / 4)) # significant bits to consider 
        nib_dx = 3 
        nibbles = [] 
        nib = 0
        for value in payload: 
            for bit_ptr in range(self.SF-1,-1, -1):
                has_ones= (value & 2**(bit_ptr)) > 0
                # print(has_ones)
                nib = nib + has_ones * 2**nib_dx     
                    
                if (nib_dx <= 0): 
                    nibbles.append(nib)
                    nib = 0
                    nib_dx = 3
                    # print("here", nib_dx)
                else: 
                    nib_dx = nib_dx - 1
                    
        if (nib_dx != 3):
            nibbles.append(nib)
        #print("Nibs:", nibbles)
        coded_nibs = self.hamming_encode(nibbles, CR) 
        
        # now select every 'SF' bits
        nsymbols = int(np.ceil((self.SF +CR) * len(payload)) / self.SF)
        symbols = []
        nib_dx = 0
        nib_offset =  (3 + CR) 
        sf_bit = self.SF - 1
        current_symbol = 0
        for nib in coded_nibs:
            for bit_ptr in range(nib_offset,-1,-1):
                has_ones= (nib & 2**(bit_ptr)) > 0
                current_symbol = current_symbol + has_ones*2**sf_bit
                
                # current_symbol = 
                if sf_bit <= 0:
                    sf_bit = self.SF - 1
                    symbols.append(current_symbol)
                    current_symbol = 0
                else: 
                    sf_bit = sf_bit - 1
                    
        # payload length may not match, check for incomplete 
        if (sf_bit != self.SF - 1 and current_symbol != 0):
            symbols.append(current_symbol)
        return symbols 
        
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
        # spsym = N*10
        
        T = spsym/Fs
        k = BW/(T) 
        f_start = -BW/2 
        pi = math.pi 
    
        # keep track of the past phase offset for phase continuity 
        ph_start = 0
    
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
            #queue = queue * np.exp(1j * ((2 * math.pi * freq_shift * self.t)+ leftover_phase - np.angle(queue[0])))  
            # out_sig = out_sig * np.exp(1j*(-np.angle(out_sig[0])+ ph_start )) # 6/17/24
            sym_out=np.append(sym_out,out_sig)
            ph_start = np.angle(out_sig[-1]) # 6/17/24
            # print(ph_start," PH")
            '''
            plt.figure(1)
            plt.specgram(out_sig) 
            plt.title('out')
            plt.show()
            '''
        # plt.figure(1)
        # plt.plot(sym_out)
        # plt.show()
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