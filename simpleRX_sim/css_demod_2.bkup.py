
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
        self.SNR_win =self.WINDOW_SIZE 
        self.PKT_SNR = []
        self.CR = CR
        self.hamming = Hamming(SF, CR)
        self.QUEUE_CNT = 0
        
    def css_demod(self, my_channel, queue, output): 
        if (len(self.noise_sig) == 0):
            self.noise_sig =  self.PREV_QUEUE[:self.SNR_win]     # prevent nan SNR measurement
        # print("Decoded leftlen:", len(self.LEFTOVER))
        print("Starting demod on a new queue", self.QUEUE_CNT)
        self.QUEUE_CNT = self.QUEUE_CNT + 1
        # print(self.PACKETS_DECODED)
        freq_shift = self.doppler_cum
        # print("Cum dop:", freq_shift)
        self.t= np.linspace(0, len(queue)/self.FS, len(queue))
        # queue = np.array(queue)
        leftover_phase = 0
        if (len(self.LEFTOVER) > 0):
            leftover_phase = np.angle(self.LEFTOVER[-1])
        # if (freq_shift != 0):
        queue = queue * np.exp(1j * ((2 * math.pi * freq_shift * self.t)+ leftover_phase - np.angle(queue[0])))   #tentatively adding leftover phase for random errors 
                
        self.PREV_QUEUE = np.concatenate((self.LEFTOVER, queue))
        PQ_TMP = self.PREV_QUEUE
        YY = len(self.LEFTOVER)
        self.LEFTOVER = []
        while (len(self.PREV_QUEUE) > self.WINDOW_SIZE): 
        # while ((len(self.PREV_QUEUE) > len(self.REF_PREAMBLE)+4*self.WINDOW_SIZE) > self.WINDOW_SIZE): 
            #print(len(self.PREV_QUEUE))
            #print("PP:", self.PREAMBLE_DEMOD == False , len(self.PREV_QUEUE) > (len(self.REF_PREAMBLE) + self.WINDOW_SIZE*2))
            if (self.PREAMBLE_DEMOD == False and len(self.PREV_QUEUE) > (len(self.REF_PREAMBLE) + self.WINDOW_SIZE*2)and self.PACKET_DETECTED):
                #print("Demod preamb at start of queue")
                # self.demod_preamble()  
                stat = self.demod_preamble()
                if stat == 1: 
                    self.PACKET_DETECTED = False              
                    # self.doppler_cum = 0
            elif self.PREAMBLE_DEMOD and self.PACKET_DETECTED:
                '''
                Decode symbols until the packet is exhausted 
                '''                
                if (self.PACKETS_DECODED < self.PACKET_LEN):
                    sym = self.symbol_demod()
                    #print(sym)
                    self.OUTPUT_PKT.append(sym)
                    self.PACKETS_DECODED = self.PACKETS_DECODED+1
                    # print("Decoded prevlen:", len(self.PREV_QUEUE))
                elif (self.END_COUNTER < len(self.END_DELIMETER)): 
                    # count the end delimeter 
                    check_sym = self.symbol_demod()
                    if (self.END_DELIMETER[self.END_COUNTER] != check_sym):
                        print(self.END_DELIMETER[self.END_COUNTER], check_sym)
                        print("Bad delimeter.")             
                        print("Dop err",self.doppler_cum)
                        '''
                        plt.figure(3)
                        plt.plot(PQ_TMP)
                        plt.axvline(x = YY, color = 'r') 
                        plt.show()
                        '''
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
            elif ((len(self.PREV_QUEUE)) >  1*(len(self.REF_PREAMBLE)+self.WINDOW_SIZE*3+  1001)):  # ((num_preamble)*self.WINDOW_SIZE)+self.WINDOW_SIZE*2 + ind +  nshifts
                # if :
                    # print("Looping.")
                    # break;
            
                possible_idx = np.array([])
                # possible_idx = self.pkt_detection(self.PREV_QUEUE, SF=self.SF, BW=20000, FS=200000, num_preamble=4)
                # print(len(self.PREV_QUEUE
                possible_idx = self.pkt_detection(self.PREV_QUEUE, SF=self.SF, BW=self.BW, FS=self.FS, num_preamble=9)
                possible_idx = list(set(possible_idx))
                possible_idx.sort()
                # do fine checks around the possible indices, only allow 1 pkt per frame 
                off = 2
                peakdx = 0
                dp = 0
                print("Length of pos dx:",len(possible_idx))
                # if (len(possible_idx) > 1):
                    # plt.figure(1)
                    # plt.plot(self.PREV_QUEUE)
                    # for p in possible_idx :
                        # plt.axvline(x = p, color = 'r')  
                    # plt.show()
                '''
                if (len(possible_idx) < 1):
                    print("Set to end.")
                    self.PREV_QUEUE = []
                    #self.PREV_QUEUE = np.array([])
                '''    
                # no longer need to do xcorr here - should be precise anyways! 
                if len(possible_idx) > 0:          
                    SNR_win = self.WINDOW_SIZE          
                    # if there are previous samples, use them to measure SNR [TO DO: this is not quite seamless]
                    if (len(self.PREV_QUEUE[:possible_idx[0]]) >= SNR_win):
                        # noise_power = np.sqrt(np.mean((self.PREV_QUEUE[(possible_idx[0]-self.WINDOW_SIZE):possible_idx[0]])**2))
                        self.noise_sig = self.PREV_QUEUE[(possible_idx[0]-SNR_win):possible_idx[0]]
                    
                    # print("PREV_QUEUE_LEN,", len(self.PREV_QUEUE),possible_idx[0])
                    self.PREV_QUEUE = self.PREV_QUEUE[possible_idx[0]:]
                    self.PACKET_DETECTED = True
                    # print("PREV_QUEUE_LEN,", len(self.PREV_QUEUE))        
                else: 
                    self.PACKET_DETECTED = False
                #print(len(self.PREV_QUEUE), (len(self.REF_PREAMBLE) + self.WINDOW_SIZE*2))
                # if self.PACKET_DETECTED and len(self.PREV_QUEUE) > (len(self.REF_PREAMBLE) + self.WINDOW_SIZE*2): 
                if self.PACKET_DETECTED and len(self.PREV_QUEUE) > (len(self.REF_PREAMBLE) + self.WINDOW_SIZE*4+1000): 
                    '''
                    Get packet count; apply doppler correction using reference chirps
                    
                    TODO: 
                    '''      
                    #print("Demod preamble within queue")
                    stat = self.demod_preamble()
                    if stat == 1: 
                        self.PACKET_DETECTED = False
                        self.PREV_QUEUE = self.PREV_QUEUE[1:] # this prevents it from stalling in a loop
                        # self.doppler_cum = 0
                else:                   
                    # self.TMP_QUEUE = self.PREV_QUEUE[len(self.PREV_QUEUE) - (len(self.REF_PREAMBLE)*4 + self.WINDOW_SIZE*2):]   #todo optimize this. this slows down our code
                    self.TMP_QUEUE = self.PREV_QUEUE
                    self.PREV_QUEUE = []
                # else:
                    # print("Not long enough.")
            else:
                break
        # if (len(self.TMP_QUEUE)) > 0:
            # self.LEFTOVER = self.TMP_QUEUE
            # self.TMP_QUEUE = []
        # else: 
            # self.LEFTOVER = self.PREV_QUEUE     
        self.LEFTOVER = self.PREV_QUEUE
        self.PREV_QUEUE = []
    
    def demod_preamble(self):
        '''
        Get dopppler from the preamble, and decode packet length. 
        '''            
        self.END_COUNTER = 0       
        # self.doppler_cum = 0 # reset the doppler shift for the entire packet (?)             
        # Doppler correction for the whole packet   
        tp = np.linspace(0, (len(self.REF_PREAMBLE)+2*self.WINDOW_SIZE)/self.FS, len(self.REF_PREAMBLE)+2*self.WINDOW_SIZE)
        
        
        # # plt.figure(1)
        # # plt.plot(self.PREV_QUEUE[:len(self.REF_PREAMBLE)])
        # # plt.show()
        # freq_shift = self.get_doppler_preamble(self.REF_PREAMBLE, self.PREV_QUEUE[:len(self.REF_PREAMBLE)],self.FD_MAX,self.FD_COARSE,tp)
        preamb_with_sync = np.concatenate((self.REF_PREAMBLE, np.conjugate(self.sym_2_css([0,0], self.N, self.SF, self.BW, self.FS))))
        preamb_candidate = self.PREV_QUEUE[:len(self.REF_PREAMBLE)+self.WINDOW_SIZE*2]
        freq_shift1 = self.get_doppler_preamble(preamb_with_sync, preamb_candidate,self.FD_MAX,self.FD_COARSE,tp)
        # self.doppler_cum = self.doppler_cum + freq_shift1
        # freq_shift = 0 #[remove this]
        
        # t = np.linspace(0, len(self.PREV_QUEUE)/self.FS, len(self.PREV_QUEUE))
        preamb_candidate = preamb_candidate * np.exp(1j * 2 * math.pi * freq_shift1 * tp)
        
        # freq_shift = self.get_doppler_preamble(self.REF_PREAMBLE, self.PREV_QUEUE[:len(self.REF_PREAMBLE)],self.FD_COARSE+2*self.FD_FINE,self.FD_FINE,tp)        
        # preamb_with_sync = np.concatenate((self.REF_PREAMBLE, np.conjugate(self.sym_2_css([0,0], self.N, self.SF, self.BW, self.FS))))
        freq_shift2 = self.get_doppler_preamble(preamb_with_sync,preamb_candidate ,self.FD_COARSE+2*self.FD_FINE,self.FD_FINE,tp)        
        
        freq_shift = freq_shift1 + freq_shift2
        
        # plt.figure(7)
        # plt.plot(self.PREV_QUEUE[:self.WINDOW_SIZE])
        # plt.plot(self.sym_2_css([1], self.N, self.SF, self.BW, self.FS))
        # plt.title('sync check')
        # plt.show()
        
        # print('Doppler preamble Frequency shift,', freq_shift)
        # check to make sure we can demodulate the preamble. If not, something went wrong with our sync and we need to try another pt 
        ts = np.linspace(0, self.WINDOW_SIZE/self.FS, self.WINDOW_SIZE)
        for symdx in range(0,7):
            # print(self.symbol_demod_sig(self.PREV_QUEUE[self.WINDOW_SIZE*(symdx):self.WINDOW_SIZE*(symdx+1)]))
            sym_win = self.PREV_QUEUE[self.WINDOW_SIZE*(symdx):self.WINDOW_SIZE*(symdx+1)] * np.exp(1j * 2 * math.pi * freq_shift * ts)
            if (self.symbol_demod_sig(sym_win) != 0):
                # print("Sync or doppler correction is incorrect.",self.symbol_demod_sig(sym_win))
                return 1
        
        
        # freq_shift = 0 #[remove this]
        t = ts = np.linspace(0, len( self.PREV_QUEUE)/self.FS, len( self.PREV_QUEUE))
        self.PREV_QUEUE = self.PREV_QUEUE * np.exp(1j * 2 * math.pi * freq_shift * t)    
        # print("PREAMB freq shift:" , freq_shift)
        self.OUTPUT_PKT = []
        self.doppler_cum = self.doppler_cum + freq_shift
        
        # self.PREV_QUEUE = self.PREV_QUEUE[len(self.REF_PREAMBLE):] 
        self.PREV_QUEUE = self.PREV_QUEUE[2*self.WINDOW_SIZE+len(self.REF_PREAMBLE):]            
        
        # SNR Measurement after doppler correction (for filtering) 
        self.noise_sig = scipy.signal.decimate(self.noise_sig, 10)
        # noise_sig = self.noise_sig
        SNR_win = self.WINDOW_SIZE 
        noise_power = np.sqrt(np.mean(np.real(self.noise_sig)*np.real(self.noise_sig)+ np.imag(self.noise_sig)*np.imag(self.noise_sig)))**2 / (SNR_win*self.UPSAMP)
        signal_sig = self.PREV_QUEUE[:SNR_win]
        # signal_power = np.sqrt(np.mean((self.PREV_QUEUE[possible_idx[0]:(possible_idx[0]+self.WINDOW_SIZE)])**2))
        signal_sig = scipy.signal.decimate(signal_sig, 10)
        # print(len(signal_sig))
        signal_power = np.sqrt(np.mean(np.real(signal_sig)*np.real(signal_sig)+np.imag(signal_sig)*np.imag(signal_sig)))**2 / (SNR_win*self.UPSAMP)
                       
        current_SNR = 20 * np.log10( (signal_power-noise_power) / noise_power) 
        self.PKT_SNR=current_SNR
        print("SNR For this packet: ", current_SNR)     
        # print(colored(("SNR For this packet: ", current_SNR), 'green', 'on_red'))
                
        pkt_len_2 = self.symbol_demod()
                       
        pkt_len_1 = self.symbol_demod()
                
        #[TODO] : What if other candidate IDXs lie within this, and this sync is just a bit off? Need to not remove the pts then
        self.PACKET_LEN =  pkt_len_2 + self.N*pkt_len_1
        if (self.PACKET_LEN > 1000):
            print(" ===== found pkt length too long ", self.PACKET_LEN)
            self.PREAMBLE_DEMOD = False 
            self.PACKET_DETECTED = False
            '''
            plt.figure(2)
            plt.plot(self.PREV_QUEUE)
            plt.show()     
            '''
        else: 
            #print("Packet length: ", self.PACKET_LEN)
            self.PACKETS_DECODED = 0      
            self.PREAMBLE_DEMOD = True
            
    def symbol_demod_sig(self,sig_tmp):
        '''
        Demodulates a CSS symbol and returns the frequency bin 
        at which the symbol appears.
        '''
        #sig_tmp = self.PREV_QUEUE[:self.WINDOW_SIZE]
                
        trans_upchirp = np.conjugate(self.UPCHIRP)
        dechirped = sig_tmp * trans_upchirp
        dechirped = np.squeeze(dechirped)
        data_fft = abs(fft.fft(dechirped)).transpose()
        
        dechirped = np.concatenate((data_fft[:int(self.N/2)+1], \
                    data_fft[int(self.N/2 + int(self.UPSAMP-1)*self.N + 1):]))
                 
        freq_bin = np.argmax(dechirped)    
        return freq_bin 
        
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
        
        freq_shift = self.get_doppler(self.sym_2_css([freq_bin], self.N, self.SF, self.BW, self.FS), self.PREV_QUEUE[:self.WINDOW_SIZE],self.FD_FINE*10,self.FD_FINE,self.tr) 
        freq_shift = 0 #[remove this]
        self.t= np.linspace(0, len(self.PREV_QUEUE)/self.FS, len(self.PREV_QUEUE))
        self.PREV_QUEUE = self.PREV_QUEUE * np.exp(1j * 2 * math.pi * freq_shift * self.t)    
        self.PREV_QUEUE = self.PREV_QUEUE[self.WINDOW_SIZE:]        
        return freq_bin    
              
    def DC_gen(self,SF, BW, Fs):
        sample_factor = int(Fs / BW)

        chirp_size = 2 ** SF
        symbol_length = sample_factor * chirp_size
        freq_shift = Fs / symbol_length
        symbol_time = 1 / freq_shift

        init_freq = -BW / 2
        final_freq = (BW / 2) - (freq_shift / sample_factor)
        t_arr = np.linspace(0, symbol_time - symbol_time/symbol_length, int(symbol_length))
        real = scipy.signal.chirp(t_arr, f0=init_freq,  f1=final_freq, t1=t_arr[-1], method='linear', phi=90)
        imag = scipy.signal.chirp(t_arr, f0=init_freq,  f1=final_freq, t1=t_arr[-1], method='linear', phi=180)

        return real + 1j * imag    
        
        
    def get_doppler_preamble_bkup(self, reference, rx_data, doppler_freq,step,t): 
        '''
        Corrects doppler shift based on the known preamble 
        '''       
        # sig_tmp = rx_data[:self.WINDOW_SIZE]
        # trans_upchirp = np.conjugate(reference[:self.WINDOW_SIZE])
        
        sig_tmp = rx_data[:len(reference)]
        trans_upchirp = np.conjugate(reference)
        
        
        dechirped = sig_tmp * trans_upchirp
        dechirped = np.squeeze(dechirped)
        
        
        data_fft = abs(fft.fft(dechirped)).transpose()
        
        
        dechirped = np.concatenate((data_fft[:int(self.N/2)+1], \
                    data_fft[int(self.N/2 + int(self.UPSAMP-1)*self.N + 1):]))
                 
        freq_bin = np.argmax(dechirped)
        # plt.figure(1)
        # plt.plot(dechirped)
        # plt.show()
        # print("frequency bin:", freq_bin)
        
        freq_bin_shift = freq_bin / (self.N*self.UPSAMP) * self.BW 
        freq_bin_shift = 0 # [remove this]
        
        # doppler needs to be updated at queue boundaries! 
        # self.doppler_cum = self.doppler_cum + freq_bin_shift
        
        return freq_bin_shift
        
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
     
    def get_doppler_bkup(self, reference, rx_data, doppler_freq,step,t): 
        '''
        Corrects doppler shift based on the known preamble 
        '''       
        sig_tmp = rx_data[:self.WINDOW_SIZE]
        trans_upchirp = np.conjugate(self.UPCHIRP)
        dechirped = sig_tmp * trans_upchirp
        dechirped = np.squeeze(dechirped)
        data_fft = abs(fft.fft(dechirped)).transpose()
        
        dechirped = np.concatenate((data_fft[:int(self.N/2)+1], \
                    data_fft[int(self.N/2 + int(self.UPSAMP-1)*self.N + 1):]))
                 
        freq_bin = np.argmax(dechirped)
        
        freq_bin_shift = freq_bin / (self.N*self.UPSAMP) * self.BW 
        freq_bin_shift = 0 # [remove this]
        
        # doppler needs to be updated at queue boundaries! 
        # self.doppler_cum = self.doppler_cum + freq_bin_shift
        return freq_bin_shift
        
    def get_doppler(self, reference, rx_data, doppler_freq,step,t): 
        '''
        Corrects doppler shift based on the known preamble 
        '''       
        f_d = np.arange(-doppler_freq, doppler_freq, step)
        xcorr_arr = []
        rx_data = rx_data/max(rx_data)
        
        trans_upchirp = np.conjugate(self.UPCHIRP)
        
        reference = reference * trans_upchirp
        reference = np.squeeze(reference)
        data_fft = abs(fft.fft(reference)).transpose()      
        reference = data_fft 
        
        rx_data = rx_data * trans_upchirp
        rx_data = np.squeeze(rx_data)
        data_fft = abs(fft.fft(rx_data)).transpose()        
        rx_data = data_fft
        
        for f in f_d: 
            t = np.linspace(0, len(rx_data)/self.FS, len(rx_data))
            data_shifted = fft.ifft(fft.fft(rx_data) * np.exp(1j * 2 * math.pi * f * t))
            xcorr_val = np.squeeze(abs(np.correlate(data_shifted, reference))) 
            
            xcorr_arr.append(xcorr_val)
                       
        freq_bin_shift = f_d[np.argmax(xcorr_arr)]
        # freq_bin_shift = 0 # [remove this]
        # doppler needs to be updated at queue boundaries! 
        # self.doppler_cum = self.doppler_cum + freq_bin_shift
        
        return freq_bin_shift
        
    def pkt_detection(self, Rx_Buffer, SF, BW, FS, num_preamble):
        upsampling_factor = int(FS / BW)
        # print("Upsampling factor:", upsampling_factor)
        N = int(2 ** SF)
        num_preamble -= 1  # need to find n-1 total chirps (later filtered by sync word)

        # Rx_Buffer_bkup = Rx_Buffer 
        # Rx_buffer = Rx_Buffer[(((num_preamble)*self.WINDOW_SIZE)+self.WINDOW_SIZE*2) : (len(Rx_Buffer) -((num_preamble)*self.WINDOW_SIZE)+self.WINDOW_SIZE*2)]

        DC_upsamp = np.conjugate(self.sym_2_css([0],  self.N, self.SF, self.BW, self.FS))
        # Preamble Detection
        ind_buff = np.array([])
        count = 0
        Pream_ind = np.array([], int)

        loop = 0
        # for off in range(5):
        # for off in range(3):
        for off in range(3):
            offset = off * upsampling_factor * N // 3
            loop = Rx_Buffer.size // (upsampling_factor * N) - 1
            for i in range(loop):
                temp_wind_fft = abs(
                    np.fft.fft(Rx_Buffer[(i * upsampling_factor * N) + offset:
                                         ((i + 1) * upsampling_factor * N) + offset] * DC_upsamp, axis=0))
                temp_wind_fft_idx = np.concatenate(
                    [np.arange(0, N // 2), np.arange(N // 2 + (upsampling_factor - 1) * N, upsampling_factor * N)])
                # temp_wind_fft_idx = np.arange(0, upsampling_factor * N)
                
                temp_wind_fft = temp_wind_fft[temp_wind_fft_idx]
                
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
        
        # Rx_Buffer = Rx_Buffer_bkup
        # Pream_ind = Pream_ind + (((num_preamble)*self.WINDOW_SIZE)+self.WINDOW_SIZE*2) 
        '''Added portion'''
        out_preamb = []
        off = 50 # number of samples to compensate for correlation offsets
        for ind in Pream_ind: 
            # possible_samples = Rx_Buffer[ind-1*self.PREAMBLE_SIZE:ind+self.PREAMBLE_SIZE*6+off]
            # possible_samples = Rx_Buffer[ind-1*self.PREAMBLE_SIZE:ind+self.PREAMBLE_SIZE*6+off]
            # [detected, peakdx] = self.checkXCORR(possible_samples)       
            
            win_size = int(len(self.REF_PREAMBLE) )
            # ind = (ind-(1*win_size)-off) + peakdx -2 #+ self.WINDOW_SIZE # NEEDS TO BE REMOVED< ignore sync word
            ind = ind 
            # out_preamb.append(ind) 
            nshifts = 1000
            shifts = np.arange(-nshifts, nshifts)
            # shifts = [0]
            # shifts = [0]
            max_shift_arr = []
            # print("check shifts")
            # if (((num_preamble)*self.WINDOW_SIZE)+self.WINDOW_SIZE*2 + ind +  nshifts < len(Rx_Buffer) and (ind- ((num_preamble)*self.WINDOW_SIZE)+self.WINDOW_SIZE*2  +  nshifts) > 0):         
            if (True):
                for shift in shifts:
                    
                    # # do sync word detection here
                    ind_win = ((num_preamble)*self.WINDOW_SIZE) + ind + shift
                    # win1 = Rx_Buffer[ind_win-self.WINDOW_SIZE:ind_win]
                    # win2 = Rx_Buffer[ind_win:ind_win+self.WINDOW_SIZE]   # transmitted downchirp 
                    # print("Current shift:",shift, ind_win)
                    # win1 = Rx_Buffer[ind_win-2*self.WINDOW_SIZE:ind_win]
                    # win1 = self.sym_2_css([0,0], self.N, self.SF, self.BW, self.FS)
                    # win2 = Rx_Buffer[ind_win:ind_win+2*self.WINDOW_SIZE]
                    # win2 = Rx_Buffer[ind_win:ind_win+2*self.WINDOW_SIZE]
                    # win1 =  Rx_Buffer[ind_win-self.WINDOW_SIZE*2:ind_win]
                    
                    win2 = np.concatenate(( np.conjugate(self.sym_2_css([0], self.N, self.SF, self.BW, self.FS)),self.sym_2_css([0], self.N, self.SF, self.BW, self.FS)))
                    # win2 = np.concatenate(( self.sym_2_css([0], self.N, self.SF, self.BW, self.FS),np.conjugate(self.sym_2_css([0], self.N, self.SF, self.BW, self.FS))))
                    win1 =  Rx_Buffer[ind_win-self.WINDOW_SIZE:ind_win+self.WINDOW_SIZE]
                    # win2 = Rx_Buffer[ind_win:ind_win+self.WINDOW_SIZE]
                    # win1 =  Rx_Buffer[ind_win-self.WINDOW_SIZE:ind_win]
                    
                    # plt.figure(1)
                    # plt.specgram(win1)
                    
                    # plt.figure(2)
                    # plt.specgram(win2)
                    # plt.show()
                    # plt.figure(3)
                    # plt.plot(Rx_Buffer[ind_win-self.WINDOW_SIZE:ind_win+self.WINDOW_SIZE*2])
                    # plt.show()
                    
                    temp_wind_fft = abs(
                            np.fft.fft(win1 * win2, axis=0))
                    
                    # plt.figure(4)
                    # plt.plot(temp_wind_fft)
                    # plt.show()
                    temp_wind_fft_idx = np.concatenate(
                        [np.arange(0, N // 2), np.arange(N // 2 + (upsampling_factor - 1) * N, upsampling_factor * N)])
                    temp_wind_fft = temp_wind_fft[temp_wind_fft_idx]
                    b = np.argmax(temp_wind_fft)                           
                    
                    # b = b + 10 
                    #print("b", b, temp_wind_fft[b])
                    max_shift_arr.append(temp_wind_fft[b])    
                    # if b == 0 : 
                        # print("b", b, temp_wind_fft[b])
                        # detected = True
                        # # max_shift_arr.append(temp_wind_fft[b])                    
                    # # else : 
                        # # max_shift_arr.append(0)
                        # detected = False 
                        # pass
                        
                # plt.figure(5)
                # plt.plot(max_shift_arr)
                # plt.show()
                            
                if (len(max_shift_arr) > 0):
                    # ind = ind + shifts[np.argmax(max_shift_arr)] + 1*self.WINDOW_SIZE#- len(self.REF_PREAMBLE)
                    ind = ind + shifts[np.argmax(max_shift_arr)] #+ 2*self.WINDOW_SIZE
                # if detected: 
                    # ind = ind + 2*self.WINDOW_SIZE # compensate for sync word offset
                    # ind = ind
                    print("PACKET DETECTED!")
                    out_preamb.append(ind) 
                    break; # tentatively adding this, may want to remove 
            else: 
                # plt.figure(1)
                # plt.plot(Rx_Buffer)
                # plt.axvline(x = ind, color = 'r') 
                # plt.show()
                print("Not enough points for packet detection ... ")
            
        Pream_ind = out_preamb        
        
        # print('Num preamble indices total:', len(Pream_ind))
        
        return Pream_ind    
        
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
        
    def create_upchirp(self):
        '''
        Create an upchirp which is used to demodulate the symbols 
        '''
        return self.sym_2_css([0], self.N, self.SF, self.BW, self.FS)
        
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