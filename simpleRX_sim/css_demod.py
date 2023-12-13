'''
css_demod :: Jackie Schellberg :: 10/19/2023 

Code to demodulate received samples from the USRP according to CSS 

my_channel : unused 

input_queue : rx complex samples 

output_queue : demodulated symbols 
'''

from matplotlib import pyplot as plt
import time
import math
import numpy as np
from scipy import signal, fft
from scipy.signal import find_peaks
from statistics import mean

import scipy.io

class CssDemod:
    def __init__(self, N, UPSAMP,PREAMBLE_SIZE,END_DELIMETER): 
        '''
        Initialize our CSS demodulator, keeping track of the state. 
        '''
        # we need to keep track if  we have detected the preamble 
        self.PACKET_DETECTED = False 
        
        # Number of samples per symbol = 2^SF 
        self.N = int(N)
        
        # Threshold envelope; at what power level do we expect to see a packet arrive? 
        # For low power scenario, this will have to be substituted 
        #self.DB_THRESH = -13 # simulation with .005 noise voltage
        #self.DB_THRESH = -25 # for simulation without any noise; should be perfect
        #self.DB_THRESH = -27
        #self.DB_THRESH = -20
        #self.DB_THRESH = -7
        #self.DB_THRESH = -30
        self.DB_THRESH = -33.4
        
        # Upsampling rate 
        self.UPSAMP = UPSAMP 
        
        # Preamble size 
        #self.PREAMBLE_SIZE = UPSAMP * N * 4
        self.PREAMBLE_SIZE = PREAMBLE_SIZE
        
        # Window size, the size of each chirp 
        self.WINDOW_SIZE = UPSAMP * N
        
        # This is a list to allocate samples as we identify them 
        self.PREV_QUEUE = []
        
        # Redundant; just makes a reference chirp 
        self.UPCHIRP = self.create_upchirp()
                
        # this will be updated during packet detection 
        self.PACKET_LEN = 59
        
        # keep track of how many packets we have decoded 
        self.PACKETS_DECODED = 0 
        
        # Leftover queue to be concatenated with the other 
        self.LEFTOVER = [] 
        
        # Doppler correction related 
        self.DOPPLER_CORR = 0 
        
        # Identifies end of the packet 
        self.END_DELIMETER = END_DELIMETER
        
        self.END_COUNTER = 0
        
        # just counting the number of packets decoded
        self.TOTAL_PKT_CNT = 0
        
    def checkXCORR(self, queue_list,possible_idx):
        #corr_offset = 50 # number of samples we want to precede the corrrelation, ensures we have a peak to detect in the presense of noise
        #corr_offset = self.PREAMBLE_SIZE
        corr_offset = 500
        print("PREV QUEUE LEN",len(self.PREV_QUEUE))
        # packet may not be centered, so we need to search all possible start of packet 
        #possible_idx = np.squeeze(np.argwhere(10*np.log10(abs(queue_list)) > self.DB_THRESH))
        
        '''
        plt.figure(1)
        xpts =range(0, len(queue_list))
        plt.plot(xpts, np.real(queue_list))
        
        #print(possible_idx)
        
        #print(possible_idx_diff)
        possible_idx_dx = np.argwhere(possible_idx_diff > self.WINDOW_SIZE)
        print(possible_idx_dx)
        plt.show()
        '''
        if (len(possible_idx) == 0) :
            self.PREV_QUEUE = [] 
            return 1
        
        print("First pos:", possible_idx[0])
        if (min(possible_idx) < corr_offset):
            possible_idx_diff = np.diff(possible_idx)  
            pdx_diff = np.atleast_1d(np.squeeze(np.argwhere(possible_idx_diff > self.PREAMBLE_SIZE)))
            #pdx_diff = np.argwhere(possible_idx_diff > self.PREAMBLE_SIZE)
            
            '''
            plt.figure(1)
            xpts =range(0, len(possible_idx_diff))
            plt.plot(xpts, (possible_idx_diff))
            plt.show()
            '''
            
            if pdx_diff.size > 0: 
                print("PDX PDX")
                possible_idx =  possible_idx[(pdx_diff[0]+1):] - corr_offset
            else:           
                print("No received signal, or the signal spills over into the previous queue item")                
                return 1    
        else: 
            possible_idx = possible_idx - corr_offset
        queue_list = queue_list[np.squeeze(possible_idx[0]):]
        
        if (len(queue_list) < self.PREAMBLE_SIZE*4): 
            return 1
        
        #trans_upchirp = np.conjugate(self.UPCHIRP)
        #trans_upchirp = np.conjugate(self.sym_to_data_ang( [1], self.N, self.UPSAMP))
        #trans_upchirp = self.sym_to_data_ang( [1], self.N, self.UPSAMP)
        #dechirped = sig_tmp * trans_upchirp
        
        xcorr_arr = [];
        #for i in range(1,self.WINDOW_SIZE*2) : 
        for i in range(1,self.PREAMBLE_SIZE*4) : 
        #for i in range(1,self.PREAMBLE_SIZE*4) : 
        #for i in range(1,len(queue_list) - 2*self.PREAMBLE_SIZE):
            # check the rest of the queue, append the remaining of the queue to the next one.                     
            window_1 = queue_list[ i : ( i + self.PREAMBLE_SIZE ) ]
            window_2 = queue_list[( i +self.PREAMBLE_SIZE) : ( i+2*self.PREAMBLE_SIZE)]
            
            fft_window1 = fft.fft(window_1)
            fft_window2 = fft.fft(window_2)
            #fft_window1 = [self.symbol_demod(window_1)]
            #fft_window2 = [self.symbol_demod(window_2)]
            
            #fft_window1 = window_1
            #fft_window2 = window_2
            # this scales the correlation between 0,1
            fft_window1 = (fft_window1 - np.mean(fft_window1)) / (np.std(fft_window1) * len(fft_window1))
            fft_window2 = (fft_window2 - np.mean(fft_window2)) / (np.std(fft_window2))
            
            xcorr_val = np.squeeze(abs(np.correlate(fft_window1, fft_window2)))
            xcorr_arr.append(xcorr_val) 
        
        #Na = 1 # for a moving average; minimize local peaks
        #xcorr_arr = np.squeeze(np.convolve(xcorr_arr, np.ones(Na)/Na, mode='valid'))
        xcorr_arr = np.array(xcorr_arr, dtype="object")
        
        # reject samples without any correlation 
        if (max(xcorr_arr) < .5):
            self.PREV_QUEUE = self.PREV_QUEUE[:possible_idx[-1]]
            return
        
        
        peaks, _ = find_peaks(xcorr_arr,height=.1) 
        #imax_peak = peaks[np.argmax(xcorr_arr[peaks])]
        imax_peak = np.argmax(xcorr_arr)
        '''
        if len(peaks) < 1:
            print("NO peaks.")
            return 
        '''
        argm = np.argmax(xcorr_arr)
        max_xcorr = max(xcorr_arr)
                      
        #peakdx = int(possible_idx[imax_peak]) + 2
        #peakdx = int(possible_idx[imax_peak]) 
        peakdx = int(possible_idx[0]+imax_peak ) - 7        # why the magic number 7?
        #peakdx = int(possible_idx[imax_peak]) -7
        #peakdx = int(possible_idx[imax_peak])
        print("Second pos:", peakdx)
        
        '''
        plt.figure(1)
        xpts =range(0, len(xcorr_arr))
        plt.plot(xpts, (xcorr_arr))
        plt.axvline(x = imax_peak, color = 'b', label = 'axvline - full height')
        
        plt.figure(2)
        xpts =range(0, len(self.PREV_QUEUE))
        plt.plot(xpts, (self.PREV_QUEUE))
        plt.axvline(x = possible_idx[0], color = 'r', label = 'axvline - full height')
        plt.axvline(x = peakdx, color = 'b', label = 'axvline - full height')
        plt.show()
        '''
        
        self.PREV_QUEUE = self.PREV_QUEUE[peakdx:]
       
        self.PACKET_DETECTED = True
                
        # get the length of the packet 
        self.PREV_QUEUE = self.PREV_QUEUE[2*self.PREAMBLE_SIZE:]
        self.PACKET_LEN = self.symbol_demod(self.PREV_QUEUE)
        print("Packet length: ", self.PACKET_LEN)
        self.PREV_QUEUE = self.PREV_QUEUE[self.WINDOW_SIZE:]
        self.PACKETS_DECODED = 0      
                
        return 0 
        
    def css_demod(self, my_channel, queue, output):    
        print("Starting Demod")
        self.PREV_QUEUE = np.concatenate((self.LEFTOVER, list(queue)))
        '''
        plt.figure(1)
        xpts =range(0, len(self.PREV_QUEUE))
        plt.plot(xpts, 10*np.log10(abs(self.PREV_QUEUE)))
        plt.show()
        '''
        
        while (len(self.PREV_QUEUE) > self.WINDOW_SIZE): 
            
            if (self.PACKET_DETECTED == False):
                possible_idx = np.squeeze(np.argwhere(10*np.log10(abs(self.PREV_QUEUE)) > self.DB_THRESH))
                #possible_idx = np.array(range(20, len(self.PREV_QUEUE)))
                if (possible_idx.size <= 1):
                    print("Not enough indices: only ", possible_idx.size)
                    print("Max Rx: ", max(10*np.log10(abs(self.PREV_QUEUE))))
                    self.PREV_QUEUE=[]
                    break
                
                ret = self.checkXCORR(self.PREV_QUEUE,possible_idx) 
                if (ret == 1): 
                    break
            elif (self.PACKET_DETECTED == True): 
                if (self.PACKET_LEN > self.PACKETS_DECODED):
                    freq_bin = self.symbol_demod(self.PREV_QUEUE) 
                    output.append(freq_bin)
                    self.PREV_QUEUE = self.PREV_QUEUE[self.WINDOW_SIZE:]
                    self.PACKETS_DECODED = self.PACKETS_DECODED + 1
                else:
                    # check that our end delimeter is there 
                    if (self.END_COUNTER < len(self.END_DELIMETER)): 
                        sym = self.END_DELIMETER[self.END_COUNTER]
                        curr_sym = self.symbol_demod(self.PREV_QUEUE)
                        self.PREV_QUEUE = self.PREV_QUEUE[self.WINDOW_SIZE:]
                        if (sym != curr_sym):
                                #print(sym, curr_sym)
                                print("Bad delimeter; sync might not be correct.")
                                self.END_COUNTER = 0
                                self.PACKET_DETECTED = False
                                break
                        self.END_COUNTER = self.END_COUNTER + 1
                    else: 
                        self.TOTAL_PKT_CNT = self.TOTAL_PKT_CNT + 1
                        print("TOTAL PACKET COUNT: ", self.TOTAL_PKT_CNT)
                        self.END_COUNTER = 0
                        self.PACKET_DETECTED = False
                        self.PACKETS_DECODED = 0
        self.LEFTOVER = self.PREV_QUEUE
        
    def get_doppler(self, up_chirp, window1): 
        '''
        Corrects doppler shift based on the known preamble 
        '''       
        dechirped_shift = window1 * np.conjugate(up_chirp) 
        dechirped_fft =fft.fft(dechirped_shift) 
        dechirped_fft_shift = np.argmax(abs(dechirped_fft))
        
        freq_1_shift = self.sym_to_data_ang( [0], self.N, self.UPSAMP) * np.conjugate(up_chirp) 
        freq_1_shift_fft = fft.fft(freq_1_shift)
        freq_1_shift_fft_shift = np.argmax(freq_1_shift_fft)
        
        freq_bin_shift = dechirped_fft_shift - (freq_1_shift_fft_shift +1)
                
        return freq_bin_shift
        
    def symbol_demod(self, rx_sig):
        '''
        Demodulates a CSS symbol and returns the frequency bin 
        at which the symbol appears.
        '''
        sig_tmp = rx_sig[:self.WINDOW_SIZE]
        trans_upchirp = np.conjugate(self.UPCHIRP)
        dechirped = sig_tmp * trans_upchirp
        dechirped = np.squeeze(dechirped)
        data_fft = abs(fft.fft(dechirped)).transpose()
        
        dechirped = np.concatenate((data_fft[:int(self.N/2)+1], \
                    data_fft[int(self.N/2 + int(self.UPSAMP-1)*self.N + 1):]))
                 
        freq_bin = np.argmax(dechirped)
        
        return freq_bin 
        
    def create_upchirp(self):
        '''
        Create an upchirp which is used to demodulate the symbols 
        '''
        return self.sym_to_data_ang([0],self.N, self.UPSAMP)
      
    def sym_to_data_ang(self, symbol,N, UPSAMP):
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