from matplotlib import pyplot as plt
import time
import math
import numpy as np
from scipy import signal, fft
from scipy.signal import find_peaks
from statistics import mean 
import pickle
import scipy.io
import scipy.signal
from Hamming import Hamming
from demod_funcs import * # all our functions for demodulating and packet detection 

class CssDemod():
    def __init__(self, N, UPSAMP,PREAMBLE_SIZE,END_DELIMETER, 
        DB_THRESH, GND_TRUTH_PKT,EXP_PAY_LEN,EXP_PKTS,SF,BW,FS,CR=0):
    #def __init__(self, N, UPSAMP,PREAMBLE_SIZE,END_DELIMETER, 
    #    DB_THRESH, GND_TRUTH_PKT=[],EXP_PAY_LEN=0,EXP_PKTS=0,SF=9,BW=20000,Fs=200000):
        '''
        Initialize our CSS demodulator, keeping track of the state. 
        '''
        self.SF = SF
        self.BW = BW
        # self.Fs = 200000
        
        # we need to keep track if  we have detected the preamble 
        self.PACKET_DETECTED = False 
        
        # Number of samples per symbol = 2^SF 
        self.N = int(N)
        
        self.DB_THRESH = DB_THRESH
        
        # Upsampling rate 
        self.UPSAMP = UPSAMP 
        # self.UPSAMP = int(FS/N)
        
        # Preamble size 
        #self.PREAMBLE_SIZE = UPSAMP * N * 4
        self.PREAMBLE_SIZE = PREAMBLE_SIZE
                        
        #self.PREAMBLE_SIZE = len(self.REF_PREAMBLE)
        
        # Window size, the size of each chirp 
        self.WINDOW_SIZE = UPSAMP * N
        
        # This is a list to allocate samples as we identify them 
        self.PREV_QUEUE = []
                               
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
        
        self.PREAMBLE_DEMOD=False 
        
        # just counting the number of packets decoded
        self.TOTAL_PKT_CNT = 0
        
        self.OUTPUT_PKT = []
        
        # maximum doppler frequency to consider
        self.FD_MAX = int(9e3)
        
        #self.FD_FINE = 100
        #self.FD_FINE = np.floor((BW/(2**SF))/2)
        #self.FD_FINE = np.floor((BW/2)//(2**SF)) 
        self.FD_FINE = .5 # this should be based off of the number of FFT points ?
        
        self.FD_COARSE = 50 # for preamble detection 
        # self.FD_COARSE = 100 # for preamble detection 
        
        # compute this based on parameters?
        # self.FS = 200000
        self.FS = FS
        # Redundant; just makes a reference chirp 
        self.UPCHIRP = self.create_upchirp()
        
        self.tr = np.linspace(0, self.WINDOW_SIZE/self.FS, self.WINDOW_SIZE)
        
        self.doppler_cum = 0
        
        self.BAD_DELIM = False
        self.TMP_QUEUE = []
        # error measuremnts
        self.GND_TRUTH_PKT=GND_TRUTH_PKT
        self.EXP_PAY_LEN=EXP_PAY_LEN
        self.EXP_PKTS=EXP_PKTS
                
        self.OUTFILE = 'tmp.pkl' 
        self.noise_sig = []     # keep track of noise power for SNR calculation
        #self.REF_PREAMBLE = self.sym_to_data_ang( [1,1,1], self.N, self.UPSAMP)
        # self.REF_PREAMBLE = self.sym_2_css([1,1,1], self.N, self.SF, self.BW, self.FS) # assuming that Tx knows this
        self.REF_PREAMBLE = self.sym_2_css([0,0,0,0,0,0,0,0], self.N, self.SF, self.BW, self.FS) # assuming that Tx knows this
        # print("Other preamble:", len(self.REF_PREAMBLE))
        #self.REF_PREAMBLE = self.REF_PREAMBLE[:int(3*len(self.REF_PREAMBLE)/4)]
        self.SNR_win =self.WINDOW_SIZE*10
        self.PKT_SNR = []
        self.CR = CR
        self.hamming = Hamming(SF, CR)
        self.QUEUE_CNT = 0
        self.possible_idx = []
        self.possible_doppler_shifts = []
        
    def css_demod(self, my_channel, queue, output): 
        if (len(self.noise_sig) == 0):
            self.noise_sig =  queue[:self.SNR_win]     # prevent nan SNR measurement
        # print("Decoded leftlen:", len(self.LEFTOVER))
        print("Starting demod on a new queue", self.QUEUE_CNT)
        self.QUEUE_CNT = self.QUEUE_CNT + 1

        freq_shift = self.doppler_cum
        # print("Cum dop:", freq_shift)
        self.t = np.linspace(0, len(queue)/self.FS, len(queue))
        # queue = np.array(queue)
        leftover_phase = 0
        if (len(self.LEFTOVER) > 0):
            leftover_phase = np.angle(self.LEFTOVER[-1])
        queue = queue * np.exp(1j * ((2 * math.pi * freq_shift * self.t)+ leftover_phase - np.angle(queue[0])))   #tentatively adding leftover phase for random errors 
                
        queue = np.concatenate((self.LEFTOVER, queue))
        self.LEFTOVER = [] # some symbols or preambles carry over into the next item 
        
        queue_idx = 0 
        # while (len(self.PREV_QUEUE) > self.WINDOW_SIZE): 
        while (queue_idx < (len(queue)-self.WINDOW_SIZE)): 
            if (self.PREAMBLE_DEMOD == False and queue_idx < (len(queue)-(len(self.REF_PREAMBLE) + self.WINDOW_SIZE*2))and self.PACKET_DETECTED):
                # demodulate the preamble, carried over from previous (idx 1)
                while (self.PREAMBLE_DEMOD == False and len(self.possible_idx) > 0):
                    possible_idx = self.possible_idx.pop(0)
                    possible_dop = self.possible_doppler_shifts.pop(0)
                    end_preamble_dx = possible_idx + len(self.REF_PREAMBLE) + self.WINDOW_SIZE*4 #todo [sub in automatic length]
                    
                    self.PREAMBLE_DEMOD, self.PACKET_LEN, freq_shift = self.demod_preamble(queue[possible_idx:end_preamble_dx], self.REF_PREAMBLE, self.WINDOW_SIZE, self.FS, -possible_dop, self.FD_COARSE, self.FD_FINE,self.N,self.SF,self.BW)
                   
                    # demod_status, packet_len, freq_shift
                    if (self.PREAMBLE_DEMOD):
                        self.END_COUNTER = 0 
                        self.PACKETS_DECODED = 0 
                        self.doppler_cum = self.doppler_cum + freq_shift
                        
                        # calculate SNR if packet detected 
                        # ds_factor = int(self.WINDOW_SIZE / (self.N * 10))  
                        ds_factor = 10
                        noise_sig = scipy.signal.decimate(self.noise_sig, ds_factor)
                        # noise_sig = self.noise_sig
                        SNR_win = self.SNR_win
                        signal_sig = queue[possible_idx:possible_idx+SNR_win]
                        signal_sig = scipy.signal.decimate(signal_sig, ds_factor)                        
                        self.PKT_SNR = self.calc_SNR(signal_sig, noise_sig)
                        print("SNR For this packet:",self.PKT_SNR)
                        
                        t_queue = np.arange(0, len(queue))/self.FS
                        queue = queue * np.exp(1j * 2 * math.pi * freq_shift * t_queue)
                if (self.PREAMBLE_DEMOD == False):
                    self.PACKET_DETECTED = False 
                else: 
                    queue_idx = possible_idx + len(self.REF_PREAMBLE) + self.WINDOW_SIZE*4
                
            elif self.PREAMBLE_DEMOD and self.PACKET_DETECTED:
                '''
                Decode symbols until the packet is exhausted 
                '''                
                if (self.PACKETS_DECODED < self.PACKET_LEN):
                    sym = self.symbol_demod_sig(queue[queue_idx:queue_idx+self.WINDOW_SIZE])
                    self.OUTPUT_PKT.append(sym)
                    self.PACKETS_DECODED = self.PACKETS_DECODED+1
                elif (self.END_COUNTER < len(self.END_DELIMETER)): 
                    # count the end delimeter 
                    check_sym = self.symbol_demod_sig(queue[queue_idx:queue_idx+self.WINDOW_SIZE])
                    if (self.END_DELIMETER[self.END_COUNTER] != check_sym):
                        print(self.END_DELIMETER[self.END_COUNTER], check_sym)
                        print("Bad delimeter.,Dop err",self.doppler_cum)      
                    self.END_COUNTER = self.END_COUNTER + 1 
                    #self.PREAMBLE_DEMOD = False
                else: 
                    self.PACKETS_DECODED = 0
                    
                    # Apply hamming decoding 
                    self.OUTPUT_PKT = self.hamming.decode(self.OUTPUT_PKT, self.CR)                    
                    
                    output = self.OUTPUT_PKT
                    #print(output)
                    self.PACKET_DETECTED = False
                    self.TOTAL_PKT_CNT  = self.TOTAL_PKT_CNT + 1
                    self.error_measurement()
                    self.PREAMBLE_DEMOD = False
                    #print("Cumulative doppler across this packet", self.doppler_cum )
                    # self.doppler_cum = 0 # no reason to propagate this if the duty cycle is too long  
                queue_idx = queue_idx + self.WINDOW_SIZE
            # elif ((queue_idx) <  3*(len(self.REF_PREAMBLE)+self.WINDOW_SIZE*3+  1001)):  # ((num_preamble)*self.WINDOW_SIZE)+self.WINDOW_SIZE*2 + ind +  nshifts
            elif ((queue_idx) <  int(self.FS/2)):   
                possible_idx = np.array([])
                possible_idx, possible_doppler_shifts = self.pkt_detection(queue[queue_idx:], SF=self.SF, BW=self.BW, FS=self.FS, num_preamble=9)
                
                # no longer need to do xcorr here - should be precise anyways! 
                if len(possible_idx) > 0:          
                    self.PACKET_DETECTED = True
                    self.PREAMBLE_DEMOD = False
                    
                    self.possible_idx = possible_idx   
                    self.possible_doppler_shifts = possible_doppler_shifts
                else: 
                    self.PACKET_DETECTED = False
                    print("No packet in this queue.")
                    # queue_idx = 0 
                break;
            else:
                print("Prev queue not big enough for pkt detection.")
                break
        
        self.LEFTOVER = queue[queue_idx:]
        self.PREV_QUEUE = []
        
    def pkt_detection(self, Rx_Buffer, SF, BW, FS,  num_preamble):
        upsampling_factor = int(FS / BW)
        # print("Upsampling factor:", upsampling_factor)
        N = int(2 ** SF)
        num_preamble -= 1  # need to find n-1 total chirps (later filtered by sync word)
        DC_upsamp = np.conjugate(self.sym_2_css([0],  N, SF, BW, FS))
        window_size = len(DC_upsamp)
        # Preamble Detection
        ind_buff = np.array([])
        count = 0
        Pream_ind = np.array([], int)

        loop = 0
        for off in range(3):
            offset = off * upsampling_factor * N // 3
            loop = Rx_Buffer.size // (upsampling_factor * N) - 1
            for i in range(loop):
                temp_wind_fft = abs(
                    np.fft.fft(Rx_Buffer[(i * upsampling_factor * N) + offset:
                                         ((i + 1) * upsampling_factor * N) + offset] * DC_upsamp, axis=0))
                # temp_wind_fft_idx = np.concatenate(
                    # [np.arange(0, N // 2), np.arange(N // 2 + (upsampling_factor - 1) * N, upsampling_factor * N)])
                # temp_wind_fft_idx = np.arange(0, upsampling_factor * N)
                
                # temp_wind_fft = temp_wind_fft[temp_wind_fft_idx]
                
                b = np.argmax(temp_wind_fft)
                if len(ind_buff) >= num_preamble:
                    ind_buff = ind_buff[-(num_preamble - 1):]
                    ind_buff = np.append(ind_buff, b)                    
                else:
                    #pass
                    ind_buff = np.append(ind_buff, b)
                    #print("Overlap")
                    
                # if ((sum(abs(np.diff(np.mod(ind_buff, N + 1)))) <= (num_preamble + 4) or
                     # sum(abs(np.diff(np.mod(ind_buff, N)))) <= (num_preamble + 4) or
                     # sum(abs(np.diff(np.mod(ind_buff, N - 1)))) <= (num_preamble + 4)) and
                        # ind_buff.size >= num_preamble - 1):    
                # if ((sum(abs(np.diff(np.mod(ind_buff, N )))) <= (num_preamble )  and
                    # ind_buff.size >= num_preamble - 1)):   
                if ((sum(abs(np.diff(np.mod(ind_buff, N + 1)))) <= (num_preamble) or
                     sum(abs(np.diff(np.mod(ind_buff, N)))) <= (num_preamble) or
                     sum(abs(np.diff(np.mod(ind_buff, N - 1)))) <= (num_preamble)) and
                        ind_buff.size >= num_preamble - 1):
                    if np.sum(np.abs(Rx_Buffer[(i * upsampling_factor * N)
                                               + offset:((i + 1) * upsampling_factor * N) + offset])) != 0:
                        count = count + 1
                        Pream_ind = np.append(Pream_ind, (i - (num_preamble - 1)) * (upsampling_factor * N) + offset)
        
        # print('Found ', count, ' Preambles', len(Pream_ind))
        if (count == 1): 
            return []
        if count >= (loop * 0.70):
            Preamble_ind = np.array([], int)
            return Preamble_ind
                
        # # Synchronization
        Pream_ind.sort()
        
        print("Coarse preamb lin:", len(Pream_ind))
        
        '''Added portion'''
        out_preamb = []
        possible_dop = []
        for ind in Pream_ind: 
            nshifts = N*int(FS/BW)
            shifts = np.arange(0, nshifts)

            max_shift_arr = []
            max_bin_arr = []
            if (((num_preamble)*window_size)+window_size*2 + ind +  nshifts < len(Rx_Buffer)  and ind >0):
                for shift in shifts:
                    
                    # # do sync word detection here
                    ind_win = ((num_preamble)*window_size) + ind + shift
                    win2 = np.concatenate(( np.conjugate(self.sym_2_css([0], N, SF, BW, FS)),self.sym_2_css([0], N, SF, BW, FS)))
                    win1 =  Rx_Buffer[ind_win-window_size:ind_win+window_size]
                    
                    temp_wind_fft = abs(
                            np.fft.fft(win1 * win2, axis=0))
                    b = np.argmax(temp_wind_fft)   
                    max_bin_arr.append(b)
                    max_shift_arr.append(temp_wind_fft[b]) 
                            
                if (len(max_shift_arr) > 0):
                    ind = ind + shifts[np.argmax(max_shift_arr)] 
                    # print("PACKET DETECTED!")
                    out_preamb.append(ind) 
                    max_bin = max_bin_arr[np.argmax(max_shift_arr)]
                    # print(max_bin)
                    
                    # for the case of negative doppler 
                    if (max_bin > ((self.N*self.UPSAMP)/2)):    
                        max_bin = -(2*(self.N*self.UPSAMP) - max_bin)
                        # print("New bin:", (self.N*self.UPSAMP),max_bin)
                    freq_shift = self.FS/2 / (self.N*self.UPSAMP) * max_bin 
                    possible_dop.append(int(freq_shift))
                    # print("Est. freq shift: ", freq_shift)
                    # break; # tentatively adding this, may want to remove 
            else: 
                print("Not enough points for packet detection ... ")
            
        Pream_ind = out_preamb        
        Dop_out = possible_dop        
        return Pream_ind, Dop_out     
 
    def demod_preamble(self, queue, ref_preamble, window_size, FS,freq_shift1,FD_COARSE,FD_FINE,N,SF,BW):
        '''
        Get dopppler from the preamble, and decode packet length. 
        '''            
        # plt.figure(1)
        # plt.specgram(queue)
        # plt.show()
        
        demod_status = False     
        
        tp = np.linspace(0, (len(ref_preamble)+2*window_size)/FS, len(ref_preamble)+2*window_size)
        
        preamb_with_sync = np.concatenate((ref_preamble, np.conjugate(self.sym_2_css([0,0], N, SF,BW,FS))))
        preamb_candidate = queue[:len(ref_preamble)+window_size*2]
        # freq_shift1 = self.get_doppler_preamble(preamb_with_sync, preamb_candidate,FD_MAX,FD_COARSE,tp)                
        # freq_shift1 = -freq_shift1
        preamb_candidate = preamb_candidate * np.exp(1j * 2 * math.pi * freq_shift1 * tp)
        # perform fine sync, on the order of bins
        freq_shift2 = self.get_doppler_preamble(preamb_with_sync,preamb_candidate, FD_COARSE+2*FD_FINE,FD_FINE,tp)        
        
        freq_shift = freq_shift1 + freq_shift2        
        
        # print('Doppler preamble Frequency shift,', freq_shift)
        # check to make sure we can demodulate the preamble. If not, something went wrong with our sync and we need to try another pt 
        ts = np.linspace(0, window_size/FS, window_size)
        for symdx in range(0,7):
            # print(self.symbol_demod_sig(self.PREV_QUEUE[self.WINDOW_SIZE*(symdx):self.WINDOW_SIZE*(symdx+1)]))
            sym_win = queue[window_size*(symdx):window_size*(symdx+1)] * np.exp(1j * 2 * math.pi * freq_shift * ts)
            if (self.symbol_demod_sig(sym_win) != 0):
                print("Sync or doppler correction is incorrect.",self.symbol_demod_sig(sym_win))
                demod_status = False
                
        t = ts = np.linspace(0, len(queue)/FS, len(queue))
        queue = queue * np.exp(1j * 2 * math.pi * freq_shift * t)   
        self.OUTPUT_PKT = []
        end_preamb_dx = 2*window_size+len(ref_preamble)
                        
        pkt_len_2 = self.symbol_demod_sig(queue[end_preamb_dx:end_preamb_dx+window_size])
                       
        pkt_len_1 = self.symbol_demod_sig(queue[end_preamb_dx+window_size:end_preamb_dx+2*window_size])
                
        #[TODO] : What if other candidate IDXs lie within this, and this sync is just a bit off? Need to not remove the pts then
        packet_len =  pkt_len_2 + N*pkt_len_1
        if (packet_len > 1000):
            print(" ===== found pkt length too long ", packet_len)
        else: 
            demod_status = True    
        # print("Packet length:", packet_len)    
        return demod_status, packet_len, freq_shift
    
    def calc_SNR(self, signal_sig, noise_sig):
        '''
        Calculates the SNR based on a measured signal with/without signal of interest.
        
        Should be downsampled to the appropriate BW before computation. 
        '''    
        
        #print("Input sigs",len(signal_sig),len(noise_sig))
        noise_power = np.sqrt(np.mean(np.real(self.noise_sig)*np.real(self.noise_sig)+ np.imag(self.noise_sig)*np.imag(self.noise_sig)))**2 / (len(signal_sig))
        
        signal_power = np.sqrt(np.mean(np.real(signal_sig)*np.real(signal_sig)+np.imag(signal_sig)*np.imag(signal_sig)))**2 / (len(signal_sig))
                       
        current_SNR = 20 * np.log10( (signal_power-noise_power) / noise_power) 
        #print(noise_power, signal_power)
        return current_SNR
        
    def symbol_demod_sig(self,sig_tmp):
        '''
        Demodulates a CSS symbol and returns the frequency bin 
        at which the symbol appears.
        '''
                
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
        return self.sym_2_css([0], self.N, self.SF, self.BW, self.FS)
        
    def sym_2_css(self, symbol, N, SF, BW, Fs): 
        '''
        sym_2_css :: Modulates a symbol between (0,N-1) to a chirp  
        
        N = 2^SF == does not represent upsampling factor (?)
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
           
        return sym_out
        
    def setErrorMeasurementFile(self,filename): 
        self.OUTFILE = filename 
                
    def error_measurement(self):
        with open(self.OUTFILE+'_'+str(self.TOTAL_PKT_CNT) + '.pkl', 'wb') as f: 
            pickle.dump([self.GND_TRUTH_PKT,  self.OUTPUT_PKT, self.TOTAL_PKT_CNT,self.PKT_SNR],f) 
    
        #print('NUM ERRORs:',sum(np.subtract(self.GND_TRUTH_PKT, self.OUTPUT_PKT)))
        print("Ground Truth:",self.GND_TRUTH_PKT)
        print("Output:      ", self.OUTPUT_PKT)
        if (len(self.OUTPUT_PKT) == len(self.GND_TRUTH_PKT)):
            print("num wrong: ",sum(np.subtract(self.OUTPUT_PKT , self.GND_TRUTH_PKT)))    
        print("TOTAL PACKET COUNT: ", self.TOTAL_PKT_CNT)
        # if (abs(sum(np.subtract(self.OUTPUT_PKT , self.GND_TRUTH_PKT))) > 0):
            # plt.figure(6)
            # plt.plot(abs(np.subtract(self.OUTPUT_PKT , self.GND_TRUTH_PKT)))
            # plt.figure(5)
            # plt.plot(self.GND_TRUTH_PKT) 
            # plt.plot(self.OUTPUT_PKT) 
            # plt.show()
            
    def get_doppler_preamble(self, reference, rx_data, doppler_freq,step,t): 
        '''
        Corrects doppler shift based on the known preamble 
        '''       
        f_d = np.arange(-doppler_freq, doppler_freq, step)
        xcorr_arr = []
        for f in f_d: 
            t = np.linspace(0, len(rx_data)/self.FS, len(rx_data))
            data_shifted = rx_data * np.exp(1j * 2 * math.pi * f * t)   
            xcorr_val = np.squeeze(abs(np.correlate(data_shifted, reference)))            
            xcorr_arr.append(xcorr_val)
                  
        # print(f_d)
        # plt.figure(1)
        # plt.plot(f_d,xcorr_arr)
        # plt.show()
        freq_bin_shift = f_d[np.argmax(xcorr_arr)]
        # freq_bin_shift = 0 # [remove this]
        
        # doppler needs to be updated at queue boundaries! 
        # self.doppler_cum = self.doppler_cum + freq_bin_shift
        
        return freq_bin_shift