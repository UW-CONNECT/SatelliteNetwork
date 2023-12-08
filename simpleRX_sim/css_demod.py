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
    def __init__(self, N, UPSAMP): 
        '''
        Initialize our CSS demodulator, keeping track of the state. 
        '''
        # we need to keep track if  we have detected the preamble 
        self.PACKET_DETECTED = False 
        
        # Number of samples per symbol = 2^SF 
        self.N = int(N)
        
        # Threshold envelope; at what power level do we expect to see a packet arrive? 
        # For low power scenario, this will have to be substituted 
        #self.DB_THRESH = -15 # simulation with .005 noise voltage
        self.DB_THRESH = -25 # for simulation without any noise; should be perfect
        #self.DB_THRESH = -27
        #self.DB_THRESH = -60
        
        # Upsampling rate 
        self.UPSAMP = UPSAMP 
        
        # Preamble size 
        self.PREAMBLE_SIZE = UPSAMP * N * 4
        
        # Window size, i.e., the size of each chirp 
        self.WINDOW_SIZE = UPSAMP * N
        
        # This is a list to allocate samples as we identify them 
        self.PREV_QUEUE = []
        
        # Redundant; just makes a reference chirp 
        self.UPCHIRP = self.create_upchirp()
        
        # What correlation peaks to we accept in noise conditions? (for xcorr)
        self.THRESH = .98
        
        # this will need to be updated during packet detection 
        #self.PACKET_LEN = 128 + 8 - 1
        self.PACKET_LEN = 59
        #self.PACKET_LEN = 50
        
        # keep track of how many packets we have decoded 
        self.PACKETS_DECODED = 0 
        
        # Leftover queue to be concatenated with the other 
        self.LEFTOVER = [] 
        
        # Doppler correction related 
        self.DOPPLER_CORR = 0 
        
        self.END_DELIMETER = [1,1,1,1,1,1,1,1]
        
    def checkXCORR(self, queue_list):
        '''
        nosigdx =  np.argwhere(abs(queue_list) > 0)
        queue_list = np.squeeze(queue_list[nosigdx])
        print(len(queue_list))
        '''
        
        '''
        plt.figure(4)
        xpts = range(0,len(self.PREV_QUEUE))
        #plt.axvline(x = peakdx, color = 'b', label = 'axvline - full height')
        plt.plot(xpts, abs(self.PREV_QUEUE))
        plt.show()
        '''
        #queue_list = queue_list[10*np.log10(abs(queue_list) )> self.DB_THRESH]
        corr_offset = 50 # number of samples we want to precede the corrrelation

        possible_idx = np.argwhere(10*np.log10(abs(queue_list)) > self.DB_THRESH)
        if (min(possible_idx) < corr_offset):
            print("No received signal, or the signal spills over into the previous queue item")
            return     
        else: 
            possible_idx = possible_idx - corr_offset
        #print(possible_idx[1:10],possible_idx.shape)
        queue_list = queue_list[np.squeeze(possible_idx[0]):]
        
        
        xcorr_arr = [];
        print((len(queue_list) - self.PREAMBLE_SIZE ))
        #for i in range(1,(len(queue_list) - self.PREAMBLE_SIZE -1)) : 
        for i in range(1,self.WINDOW_SIZE*2) : 
            #print(i)
            # check the rest of the queue, append the remaining of the queue to the next one.                     
            window_1 = queue_list[ i : ( i + self.PREAMBLE_SIZE ) ]
            window_2 = queue_list[( i +self.PREAMBLE_SIZE) : ( i+2*self.PREAMBLE_SIZE)]
            
            fft_window1 = fft.fft(window_1)
            fft_window2 = fft.fft(window_2)
                                
            fft_window1 = (fft_window1 - np.mean(fft_window1)) / (np.std(fft_window1) * len(fft_window1))
            fft_window2 = (fft_window2 - np.mean(fft_window2)) / (np.std(fft_window2))
            
            #xcorr_val = abs(np.correlate(fft_window1, fft_window2))
            #print("LENGTHs", len(window_1), len(window_2))
            #print(len(fft_window1), len(fft_window2))
            xcorr_val = np.squeeze(abs(np.correlate(fft_window1, fft_window2)))
            #print("CORR SZ:", (np.correlate(fft_window1, fft_window2)).shape)
            xcorr_arr.append(xcorr_val) 
        
        Na = 1 # for a moving average; minimize local peaks
        xcorr_arr = np.squeeze(np.convolve(xcorr_arr, np.ones(Na)/Na, mode='valid'))
        #xcorr_arr = xcorr_arr**2
        #print('XCORR: ',xcorr_arr.shape)
        peaks, _ = find_peaks(xcorr_arr,height=.1) 
        imax_peak = peaks[np.argmax(xcorr_arr[peaks])]
        
        if len(peaks) < 1:
            print("NO peaks.")
            return 
        print(len(xcorr_arr))
        argm = np.argmax(xcorr_arr)
        max_xcorr = max(xcorr_arr)
                      
        #
        print("FULL PKT LENGTH:", len(np.argwhere(abs(queue_list) > 0)))
        peakdx = int(possible_idx[imax_peak]) + 2
        print(imax_peak)
        '''
        plt.figure(4)
        xpts = range(0,len(self.PREV_QUEUE))
        plt.axvline(x = peakdx, color = 'b', label = 'axvline - full height')
        plt.plot(xpts, self.PREV_QUEUE)
        '''
        
        self.PREV_QUEUE = self.PREV_QUEUE[peakdx:]
                     
        plt.figure(3)
        xpts = range(0,len(self.PREV_QUEUE))
        #plt.axvline(x = peakdx, color = 'b', label = 'axvline - full height')
        plt.plot(xpts, self.PREV_QUEUE)
      
        plt.figure(1)
        xpts = range(0,len(xcorr_arr))
        plt.plot(xpts, xcorr_arr)
        plt.axvline(x = imax_peak, color = 'b', label = 'axvline - full height')
        plt.axvline(x = (imax_peak+2*self.PREAMBLE_SIZE), color = 'b', label = 'axvline - full height')
        xpts = range(0,len(queue_list))
        plt.plot(xpts, queue_list)
        plt.legend(['xcorr','xline', 'queue list'])
        
        
        
        ''' Generate example preamble to compare against '''
        ex_preamble = self.sym_to_data_ang([1, 1, 1, 1, 1, 1, 1, 1],self.N, self.UPSAMP)
        cmp_ex_preamble = abs(ex_preamble)/max(abs(ex_preamble))
        cmp_prev_queue =abs(self.PREV_QUEUE[:len(ex_preamble)])/max(abs(self.PREV_QUEUE[:len(ex_preamble)]))
        plt.figure(4)
        xpts = range(0,len(ex_preamble))
        plt.plot(xpts, cmp_ex_preamble) 
        plt.plot(xpts, cmp_prev_queue)
        plt.plot(xpts, cmp_ex_preamble - cmp_prev_queue)
        plt.legend(['True Preamble','Cut'])
        plt.show()
        #max_idx, max_val = max(xcorr_vals)
        # if two windows are close to identical, we have identified the preamble ! 
        #if (max_xcorr > self.THRESH):
        self.PACKET_DETECTED = True
        queue = []
        
        # get the length of the packet 
        prev_len =len(self.PREV_QUEUE)
        self.PREV_QUEUE = self.PREV_QUEUE[2*self.PREAMBLE_SIZE:]
        self.PACKET_LEN = self.symbol_demod(self.PREV_QUEUE)
        self.PREV_QUEUE = self.PREV_QUEUE[self.WINDOW_SIZE:]
        self.PACKETS_DECODED = 0
        print(prev_len - len(self.PREV_QUEUE))
        print(self.PACKET_LEN)
        
        
        ## DO DOPPLER CORRECTION HERE!
        #self.DOPPLER_CORR = self.get_doppler(self.UPCHIRP, window_1)
        #print(self.DOPPLER_CORR)
        return 
        
    def css_demod(self, my_channel, queue, output):    
        print("Starting Demod")
        self.PREV_QUEUE = np.concatenate((self.LEFTOVER, list(queue)))
        if (self.PACKET_DETECTED == False):
            self.checkXCORR(self.PREV_QUEUE)
                                                
        if (self.PACKET_DETECTED):
            #print("Packet detected")
            # [TODO] -- symbol boundaires may not be preserved across queue items; keep track of leftover samples for sync
            # now process the remaining samples and print the result )
            while (len(self.PREV_QUEUE) >= self.WINDOW_SIZE and self.PACKET_DETECTED ):
                if (self.PACKET_LEN > self.PACKETS_DECODED):
                    freq_bin = self.symbol_demod(self.PREV_QUEUE)  # to match the Matlab bin implementation
                    #output.append(freq_bin - self.DOPPLER_CORR)
                    output.append(freq_bin)
                    prev_len =len(self.PREV_QUEUE)
                    self.PREV_QUEUE = self.PREV_QUEUE[self.WINDOW_SIZE:]
                    #print(prev_len - len(self.PREV_QUEUE))
                    self.PACKETS_DECODED = self.PACKETS_DECODED + 1
                else:
                    # check that our end delimeter is there 
                    '''
                    print(self.PACKET_LEN, self.PACKETS_DECODED)
                    for sym in self.END_DELIMETER:
                        curr_sym = self.symbol_demod(self.PREV_QUEUE)
                        self.PREV_QUEUE = self.PREV_QUEUE[self.WINDOW_SIZE:]
                        if (sym != curr_sym):
                            #print(sym, curr_sym)
                            print("Bad delimeter; sync might not be correct.")
                            break
                    '''
                    self.PACKET_DETECTED = False
                    self.PACKETS_DECODED = 0
                    
            # packet may carry over into the next queue item; so keep track of it 
            if (self.PACKET_DETECTED):        
                self.LEFTOVER = self.PREV_QUEUE
            else:
                self.LEFTOVER = []
                       
        # TODO : we want to minimize this elapsed time 
        #print(time.time() - t)
        
    def get_doppler(self, up_chirp, window1): 
        '''
        Corrects doppler shift based on the known preamble 
        '''       
        dechirped_shift = window1 * np.conjugate(up_chirp) 
        dechirped_fft =fft.fft(dechirped_shift) 
        dechirped_fft_shift = np.argmax(abs(dechirped_fft))
        #print(dechirped_fft_shift)
        
        freq_1_shift = self.sym_to_data_ang( [0], self.N, self.UPSAMP) * np.conjugate(up_chirp) 
        freq_1_shift_fft = fft.fft(freq_1_shift)
        freq_1_shift_fft_shift = np.argmax(freq_1_shift_fft)
        #print(freq_1_shift_fft_shift)
        
        freq_bin_shift = dechirped_fft_shift - (freq_1_shift_fft_shift +1)
        
        #print(freq_bin_shift)
        
        return freq_bin_shift
        
    def symbol_demod(self, rx_sig):
        '''
        Demodulates a CSS symbol and returns the frequency bin 
        at which the symbol appears.
        '''
        sig_tmp = rx_sig[:self.WINDOW_SIZE]
        
        '''
        plt.figure(2)
        xpts = range(0,len(sig_tmp))
        #plt.axvline(x = peakdx, color = 'b', label = 'axvline - full height')
        plt.plot(xpts, sig_tmp)
        
        plt.figure(3)
        f, t, Sxx = signal.spectrogram(sig_tmp, 1)
        plt.pcolormesh(t, f, Sxx, shading='gouraud')
        plt.ylabel('Frequency [Hz]')
        plt.xlabel('Time [sec]')
        plt.show()
        '''
        trans_upchirp = np.conjugate(self.UPCHIRP)
        dechirped = sig_tmp * trans_upchirp
        dechirped = np.squeeze(dechirped)
        data_fft = abs(fft.fft(dechirped)).transpose()
        
        dechirped = np.concatenate((data_fft[:int(self.N/2)+1], \
                    data_fft[int(self.N/2 + int(self.UPSAMP-1)*self.N + 1):]))
                 
        freq_bin = np.argmax(dechirped)
        #freq_bin = np.argmax(data_fft)
        
        return freq_bin 
        
    def create_upchirp(self):
        '''
        Create an upchirp which is used to demodulate the symbols 
        '''
        #return self.sym_to_data_ang([-1],self.N, self.UPSAMP)
        return self.sym_to_data_ang([0],self.N, self.UPSAMP)
    
    '''
    def sym_to_data_ang(self, symbol,N, UPSAMP):
        #        Via https://github.com/mananmishra11/open-lora/blob/main/std_lora/sym_to_data_ang.py
        
        data = []
        accumulator = 0
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
                
        return data
    '''
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