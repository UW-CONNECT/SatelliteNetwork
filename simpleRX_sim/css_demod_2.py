
from matplotlib import pyplot as plt
import time
import math
import numpy as np
from scipy import signal, fft
from scipy.signal import find_peaks
from statistics import mean 
import pickle
import scipy.io

class CssDemod():
    def __init__(self, N, UPSAMP,PREAMBLE_SIZE,END_DELIMETER, DB_THRESH, GND_TRUTH_PKT=[],EXP_PAY_LEN=0,EXP_PKTS=0):
        '''
        Initialize our CSS demodulator, keeping track of the state. 
        '''
        # we need to keep track if  we have detected the preamble 
        self.PACKET_DETECTED = False 
        
        # Number of samples per symbol = 2^SF 
        self.N = int(N)
        
        self.DB_THRESH = DB_THRESH
        
        # Upsampling rate 
        self.UPSAMP = UPSAMP 
        
        # Preamble size 
        #self.PREAMBLE_SIZE = UPSAMP * N * 4
        self.PREAMBLE_SIZE = PREAMBLE_SIZE
        
        self.REF_PREAMBLE = self.sym_to_data_ang( [1,1], self.N, self.UPSAMP)
        
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
        self.LEFTOVER = np.array([])
        
        # Doppler correction related 
        self.DOPPLER_CORR = 0 
        
        # Identifies end of the packet 
        self.END_DELIMETER = END_DELIMETER
        
        self.END_COUNTER = 0
        
        # just counting the number of packets decoded
        self.TOTAL_PKT_CNT = 0
        
        self.OUTPUT_PKT = []
        
        # maximum doppler frequency to consider
        self.FD_MAX = int(9e3)
        
        #self.FD_FINE = 100
        self.FD_FINE = 2
        
        # compute this based on parameters?
        self.FS = 200000
        
        self.tr = np.linspace(0, self.WINDOW_SIZE/self.FS, self.WINDOW_SIZE)
        
        self.doppler_cum = 0
        
        self.BAD_DELIM = False
        
        # error measuremnts
        self.GND_TRUTH_PKT=GND_TRUTH_PKT.astype(int)
        self.EXP_PAY_LEN=EXP_PAY_LEN
        self.EXP_PKTS=EXP_PKTS
                
        self.OUTFILE = 'tmp.pkl'
        
    def css_demod(self, my_channel, queue, output):   
        print("Starting demod on a new queue")
        freq_shift = self.doppler_cum
        self.t= np.linspace(0, len(queue)/self.FS, len(queue))
        queue = queue * np.exp(1j * 2 * math.pi * freq_shift * self.t) 
        self.PREV_QUEUE = np.squeeze(np.concatenate((self.LEFTOVER, queue)))
        while (len(self.PREV_QUEUE) > self.WINDOW_SIZE): 
            if self.PACKET_DETECTED:
                '''
                Decode symbols until the packet is exhausted 
                '''
                if (self.PACKETS_DECODED < self.PACKET_LEN):
                    sym = self.symbol_demod()
                    #print(sym)
                    self.OUTPUT_PKT.append(sym)
                    self.PACKETS_DECODED = self.PACKETS_DECODED+1
                elif (self.END_COUNTER < len(self.END_DELIMETER)): 
                    # count the end delimeter 
                    check_sym = self.symbol_demod()
                    if (self.END_DELIMETER[self.END_COUNTER] != check_sym):
                        print(self.END_DELIMETER[self.END_COUNTER], check_sym)
                        print("Bad delimeter.")                        
                    self.END_COUNTER = self.END_COUNTER + 1 
                else: 
                    self.PACKETS_DECODED = 0
                    output = self.OUTPUT_PKT
                    #print(output)
                    self.PACKET_DETECTED = False
                    self.TOTAL_PKT_CNT  = self.TOTAL_PKT_CNT + 1
                    self.error_measurement()
                    #print("Dop:", self.doppler_cum)
                pass
            else:            
                truth_val = (10*np.log10(abs(self.PREV_QUEUE)) > self.DB_THRESH) 
                truth_len = truth_val.sum() 
                #pkt_dx = np.array(0)
                
                if truth_len == 1:
                    pkt_dx = np.argwhere(truth_val)[0]
                else:
                    pkt_dx =  np.squeeze(np.argwhere(truth_val))
                
                #pkt_dx = [np.squeeze(np.nonzero(10*np.log10(abs(self.PREV_QUEUE)) > self.DB_THRESH))]
                #pkt_dx = np.argwhere(10*np.log10(abs(self.PREV_QUEUE)) > self.DB_THRESH)
                '''     
                plt.figure(1)
                plt.plot(10*np.log10(abs(self.PREV_QUEUE)))
                plt.show()
                '''
                #print(len(10*np.log10(abs(queue))),pkt_dx.shape, len(queue[pkt_dx]))
                print(truth_len)
                #if pkt_dx.size == 0: 
                #if (len(pkt_dx) == 0):
                if truth_len < 1:                
                    self.PREV_QUEUE = []
                elif (not pkt_dx[0]+self.WINDOW_SIZE*4 < len(self.PREV_QUEUE)):
                    #self.LEFTOVER = self.PREV_QUEUE 
                    #self.PREV_QUEUE = []
                    break;
                else: 
                    '''
                    print((not self.PACKET_DETECTED))
                    print(truth_len > 0)
                    print(pkt_dx[0]+self.WINDOW_SIZE*4 < len(self.PREV_QUEUE))
                    print(pkt_dx[0])
                    print(truth_len)
                    '''
                    print("CC",pkt_dx[0]+self.WINDOW_SIZE*4 < len(self.PREV_QUEUE))
                    while((not self.PACKET_DETECTED) and len(pkt_dx) > 0 and pkt_dx[0]+self.WINDOW_SIZE*4 < len(self.PREV_QUEUE)): 
                        #[self.PACKET_DETECTED, peakdx] = self.checkXCORR(queue[pkt_dx[1]:pkt_dx[1]+self.PREAMBLE_SIZE*4])
                        corr_offset = 50
                        [self.PACKET_DETECTED, peakdx] = self.checkXCORR(self.PREV_QUEUE[pkt_dx[0]-corr_offset:(pkt_dx[0]+4*self.PREAMBLE_SIZE+corr_offset)])
                        print("Checking XCORR")
                        if self.PACKET_DETECTED: 
                            '''
                            Get packet count; apply doppler correction using reference chirps
                            
                            TODO: 
                            '''         
                            self.END_COUNTER = 0
                            self.PREV_QUEUE = self.PREV_QUEUE[(pkt_dx[0] -corr_offset+ peakdx):]
                            
                            # Doppler correction for the whole packet   
                            t = np.linspace(0, len(self.REF_PREAMBLE)/self.FS, len(self.REF_PREAMBLE))
                            freq_shift = self.get_doppler(self.REF_PREAMBLE, self.PREV_QUEUE[:self.PREAMBLE_SIZE*2],self.FD_MAX,2,t)
                            print(freq_shift)
                            t = np.linspace(0, len(self.PREV_QUEUE)/self.FS, len(self.PREV_QUEUE))
                            self.PREV_QUEUE = self.PREV_QUEUE * np.exp(1j * 2 * math.pi * freq_shift * t)    
                            
                            self.OUTPUT_PKT = []
                        
                            self.PREV_QUEUE = self.PREV_QUEUE[2*self.PREAMBLE_SIZE:]
                            pkt_len_2 = self.symbol_demod()
                                           
                            pkt_len_1 = self.symbol_demod()
                                    
                            self.PACKET_LEN =  pkt_len_2 + self.N*pkt_len_1
                            
                            print("Packet length: ", self.PACKET_LEN)
                            self.PACKETS_DECODED = 0                                  
                        else:                         
                            #print(len(self.PREV_QUEUE))
                            #self.PREV_QUEUE = self.PREV_QUEUE[pkt_dx[0]+1:]
                            #print(len(self.PREV_QUEUE))
                            if (len(pkt_dx) == 1):
                                print("Set to end.")
                                self.PREV_QUEUE = self.PREV_QUEUE[pkt_dx[0]+1:]
                            
                            pkt_dx=pkt_dx[1:]
                            #print(len(pkt_dx))
                            
                            #self.PREV_QUEUE = self.PREV_QUEUE[pkt_dx[1:]]
                            #np.delete(pkt_dx,0)
                        
        self.LEFTOVER = self.PREV_QUEUE       
            #print("Peak dx: ", peakdx)
        
    def symbol_demod(self):
        '''
        Demodulates a CSS symbol and returns the frequency bin 
        at which the symbol appears.
        '''
        sig_tmp = self.PREV_QUEUE[:self.WINDOW_SIZE]
                
        trans_upchirp = np.conjugate(self.UPCHIRP)
        dechirped = sig_tmp * trans_upchirp
        dechirped = np.squeeze(dechirped)
        data_fft = abs(fft.fft(dechirped)).transpose()
        
        dechirped = np.concatenate((data_fft[:int(self.N/2)+1], \
                    data_fft[int(self.N/2 + int(self.UPSAMP-1)*self.N + 1):]))
                 
        freq_bin = np.argmax(dechirped)
        #freq_shift = self.get_doppler(self.sym_to_data_ang([freq_bin],self.N, self.UPSAMP), self.PREV_QUEUE[:self.#WINDOW_SIZE],20,2,self.tr)
        freq_shift = self.get_doppler(self.sym_to_data_ang([freq_bin],self.N, self.UPSAMP), self.PREV_QUEUE[:self.WINDOW_SIZE],20,2,self.tr)
        self.t= np.linspace(0, len(self.PREV_QUEUE)/self.FS, len(self.PREV_QUEUE))
        self.PREV_QUEUE = self.PREV_QUEUE * np.exp(1j * 2 * math.pi * freq_shift * self.t)    
        self.PREV_QUEUE = self.PREV_QUEUE[self.WINDOW_SIZE:]        
        return freq_bin 
        
    def checkXCORR(self, possible_samples):               
        detected = False 
        peakdx = 0
        
        xcorr_arr = [];
        for i in range(1,len(possible_samples) - 4*self.PREAMBLE_SIZE) :                     
            window_1 = possible_samples[ i : ( i + self.PREAMBLE_SIZE ) ]
            window_2 = possible_samples[( i +self.PREAMBLE_SIZE) : ( i+2*self.PREAMBLE_SIZE)]
            
            # this scales the correlation between 0,1
            fft_window1 = (window_1 - np.mean(window_1)) / (np.std(window_1) * len(window_1))
            fft_window2 = (window_2 - np.mean(window_2)) / (np.std(window_2))
            
            xcorr_val = np.squeeze(abs(np.correlate(fft_window1, fft_window2)))
            xcorr_arr.append(xcorr_val) 
        
        xcorr_arr = np.array(xcorr_arr, dtype="object")
        
        '''
        plt.figure(2)
        plt.plot(xcorr_arr)
        plt.show()
        '''
        
        if len(xcorr_arr) == 0:
            return detected, peakdx
        
        # reject samples without any correlation 
        if (max(xcorr_arr) > .5):
            detected = True
            print("Packet detected.")
                
        imax_peak = np.argmax(xcorr_arr)
        
        argm = np.argmax(xcorr_arr)
        max_xcorr = max(xcorr_arr)
                      
        #peakdx = int(possible_idx[0]+imax_peak ) - 7        # why the magic number 7?
        #peakdx = int(possible_idx[0]+imax_peak )  
        peakdx = argm        
        return detected, peakdx        
        
    def get_doppler(self, reference, rx_data, doppler_freq,step,t): 
        '''
        Corrects doppler shift based on the known preamble 
        '''       
        f_d = np.arange(-doppler_freq, doppler_freq, step)
        xcorr_arr = []
        for f in f_d: 
            t = np.linspace(0, len(reference)/self.FS, len(reference))
            data_shifted = rx_data * np.exp(1j * 2 * math.pi * f * t)        
               
            xcorr_val = np.squeeze(abs(np.correlate(data_shifted, reference)))   
            xcorr_arr.append(xcorr_val)
                       
        freq_bin_shift = f_d[np.argmax(xcorr_arr)]
        
        '''
        if (freq_bin_shift > 4):
            plt.figure(1)
            xpts =range(0, len(xcorr_arr))
            plt.plot(f_d, (xcorr_arr))
            plt.axvline(x = f_d[np.argmax(xcorr_arr)], color = 'b')
            plt.title('Freq. correlation')
            plt.show()
        '''
        # doppler needs to be updated at queue boundaries! 
        self.doppler_cum = self.doppler_cum + freq_bin_shift
        
        return freq_bin_shift
        
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
    
    def create_upchirp(self):
        '''
        Create an upchirp which is used to demodulate the symbols 
        '''
        return self.sym_to_data_ang([0],self.N, self.UPSAMP)
    
    def setErrorMeasurementFile(self,filename): 
        self.OUTFILE = filename 
        
    def error_measurement(self):
        with open(self.OUTFILE+'_'+str(self.TOTAL_PKT_CNT) + '.pkl', 'wb') as f: 
            pickle.dump([self.GND_TRUTH_PKT,  self.OUTPUT_PKT, self.TOTAL_PKT_CNT],f)
    
        print('NUM ERRORs:',sum(np.subtract(self.GND_TRUTH_PKT, self.OUTPUT_PKT)))
        print("Ground Truth:",self.GND_TRUTH_PKT)
        print("Output:", self.OUTPUT_PKT)
        print("TOTAL PACKET COUNT: ", self.TOTAL_PKT_CNT)