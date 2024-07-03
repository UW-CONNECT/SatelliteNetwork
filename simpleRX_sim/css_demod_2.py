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
# from dnsamp_buff import * # fine frequency correction related
from FFO_corr import * 
class CssDemod():
    def __init__(self, N, UPSAMP,PREAMBLE_SIZE,END_DELIMETER, 
        GND_TRUTH_PKT,EXP_PAY_LEN,EXP_PKTS,SF,BW,FS,CR=0):
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
        
        # self.DB_THRESH = DB_THRESH
        
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
        self.FD_FINE = self.BW/(4*N) # this should be based off of the number of FFT points ?
        
        self.FD_COARSE = 50 # for preamble detection 
        
        # compute this based on parameters?
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
        self.REF_PREAMBLE = self.sym_2_css([0,0,0,0,0,0,0,0], self.N, self.SF, self.BW, self.FS) # assuming that Tx knows this

        self.SNR_win =self.WINDOW_SIZE*10
        self.PKT_SNR = []
        self.CR = CR
        self.hamming = Hamming(SF, CR)
        self.QUEUE_CNT = 0
        self.possible_idx = []
        self.possible_doppler_shifts = []
        
        self.POSIX_START = 0
        self.CURR_POSIX_TIME = 0 
        
        # SuperDC for activity detection 
        self.superDChigh = self.genSuperDC(12, self.BW, self.FS)
        self.superDClow = self.genSuperDC(9, self.BW, self.FS)
        
    def css_demod(self, my_channel, queue, output, queue_time): 
        print("Starting demod.")
        if (len(self.noise_sig) == 0):
            self.noise_sig =  queue[:self.SNR_win]     # prevent nan SNR measurement
        if (self.QUEUE_CNT == 0):
            self.POSIX_START = queue_time 
        self.QUEUE_CNT = self.QUEUE_CNT + 1

        freq_shift = self.doppler_cum
        self.t = np.linspace(0, len(queue)/self.FS, len(queue))
        queue = np.concatenate((self.LEFTOVER, queue))
        self.LEFTOVER = [] # some symbols or preambles carry over into the next item 
        self.CURR_POSIX_TIME = queue_time
        
        queue_idx = 0 
        # while (len(self.PREV_QUEUE) > self.WINDOW_SIZE): 
        while (queue_idx < (len(queue)-self.WINDOW_SIZE)): 
            # if (self.PREAMBLE_DEMOD == False and queue_idx < (len(queue)-(len(self.REF_PREAMBLE) + self.WINDOW_SIZE*2))and self.PACKET_DETECTED):
            if (self.PREAMBLE_DEMOD == False and queue_idx < (len(queue)-(len(self.REF_PREAMBLE) + self.WINDOW_SIZE*4))and self.PACKET_DETECTED
                and max(self.possible_idx) < (len(queue)-(len(self.REF_PREAMBLE) + self.WINDOW_SIZE*10))):
                # demodulate the preamble, carried over from previous (idx 1)
                while (self.PREAMBLE_DEMOD == False and len(self.possible_idx) > 0):
                    possible_idx = self.possible_idx.pop(0)
                    possible_dop = self.possible_doppler_shifts.pop(0)
                    
                    
                    possible_idx, possible_dop = self.check_doppler_preamble(self.N, self.FS, self.SF, self.BW, 9, self.WINDOW_SIZE, queue,possible_idx)
                    # possible_idx = possible_idx + shift 
                    # plt.figure(66)
                    # plt.plot(queue)
                    # plt.axvline(possible_idx)
                    # plt.show()
                    print("Out of doppler preamble.", len(self.possible_idx), "PDS")
                    end_preamble_dx = possible_idx + len(self.REF_PREAMBLE) + self.WINDOW_SIZE*4 #todo [sub in automatic length]
                    
                    if (end_preamble_dx > len(queue) or possible_idx < 0): # 6/5/2024 change this possible 
                        print("Breaking")
                        possible_idx = queue_idx # may need to re-do pkt detection?
                        break;
                    
                    # plt.figure(18)
                    # plt.plot(queue)
                    # plt.axvline(possible_idx)
                    # plt.show()

                    
                    # print("Attempting preamble demod. Possible idx:", possible_idx, "possible doppler:", possible_dop)
                    self.PREAMBLE_DEMOD, self.PACKET_LEN, freq_shift,idx_shift = self.demod_preamble(queue[possible_idx-3:end_preamble_dx+14], self.REF_PREAMBLE, self.WINDOW_SIZE, self.FS, -possible_dop, self.FD_COARSE, self.FD_FINE,self.N,self.SF,self.BW)
                    # self.PREAMBLE_DEMOD, self.PACKET_LEN, freq_shift,idx_shift = self.demod_preamble(queue[possible_idx:end_preamble_dx], self.REF_PREAMBLE, self.WINDOW_SIZE, self.FS, -possible_dop, self.FD_COARSE, self.FD_FINE,self.N,self.SF,self.BW)          # [TODO] 7/1/2024 remove
                    # demod_status, packet_len, freq_shift
                    if (self.PREAMBLE_DEMOD):
                        self.END_COUNTER = 0 
                        self.PACKETS_DECODED = 0 
                        # self.doppler_cum = self.doppler_cum + freq_shift
                        self.doppler_cum = freq_shift
                        
                        # calculate SNR if packet detected 
                        # ds_factor = int(self.WINDOW_SIZE / (self.N * 10))  
                        ds_factor = self.UPSAMP
                        print("Noise sig len:", len(self.noise_sig))
                        # noise_sig = scipy.signal.decimate(self.noise_sig, ds_factor)
                        noise_sig = self.noise_sig[::ds_factor]
                        
                        print("Noise sig len:", len(noise_sig))
                        # noise_sig = self.noise_sig
                        SNR_win = self.SNR_win
                        # signal_sig = queue[possible_idx:possible_idx+SNR_win]
                        tt = np.arange(0, SNR_win)/self.FS
                        signal_sig = queue[possible_idx:possible_idx+SNR_win]* np.exp(1j * 2 * math.pi * self.doppler_cum * tt).transpose()
                        
                        # signal_sig = scipy.signal.decimate(signal_sig, ds_factor)
                        signal_sig = signal_sig[::ds_factor]                        

                        self.PKT_SNR = self.calc_SNR(signal_sig, noise_sig)
                        print("SNR For this packet:",self.PKT_SNR, "Init freq shift: ", freq_shift, "Packet len:", self.PACKET_LEN)
                        
                        t_queue = np.arange(0, len(queue))/self.FS
                        # queue = queue * np.exp(1j * 2 * math.pi * freq_shift * t_queue)
                # print("Out of pkt detection")
                if (self.PREAMBLE_DEMOD == False):
                    self.PACKET_DETECTED = False
                    print("Indices incorrect.", possible_idx)
                    # queue_idx = possible_idx + len(self.REF_PREAMBLE) + self.WINDOW_SIZE*4 
                    break
                else: 
                    queue_idx = possible_idx + len(self.REF_PREAMBLE) + self.WINDOW_SIZE*4 +idx_shift-3
            elif self.PREAMBLE_DEMOD and self.PACKET_DETECTED:
                '''
                Decode symbols until the packet is exhausted 
                '''                
                if (self.PACKETS_DECODED < self.PACKET_LEN):
                    # sym = self.symbol_demod_sig(queue[queue_idx:queue_idx+self.WINDOW_SIZE])
                    tt = np.arange(0, self.WINDOW_SIZE)/self.FS
                    
                    sym = self.symbol_demod_sig(queue[queue_idx:queue_idx+self.WINDOW_SIZE]* np.exp(1j * 2 * math.pi * self.doppler_cum * tt))

                    # get a new fine-grained doppler estimate  sym_2_css(self, symbol, N, SF, BW, Fs)
                    reference = self.sym_2_css([sym], self.N, self.SF, self.BW, self.FS)
                    t_queue = np.arange(0, len(queue))/self.FS
                    freq_lin = self.BW / (self.N)
                    f_shift_new = FFO_corr_sym(queue[queue_idx:queue_idx+self.WINDOW_SIZE]* np.exp(1j * 2 * math.pi * self.doppler_cum * tt), self.SF, self.BW, self.FS, self.N, sym)
                    # f_shift_new=0 #[TODO] Remove.
                    self.doppler_cum = self.doppler_cum + f_shift_new
                    
                    self.OUTPUT_PKT.append(sym)
                    self.PACKETS_DECODED = self.PACKETS_DECODED+1

                elif (self.END_COUNTER < len(self.END_DELIMETER)): 
                    # count the end delimeter 
                    # check_sym = self.symbol_demod_sig(queue[queue_idx:queue_idx+self.WINDOW_SIZE])
                    tt = np.arange(0, self.WINDOW_SIZE)/self.FS
                    check_sym = self.symbol_demod_sig(queue[queue_idx:queue_idx+self.WINDOW_SIZE]* np.exp(1j * 2 * math.pi * self.doppler_cum * tt))
                                        
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
            elif ((queue_idx) <  int(self.FS/2)):   # this should be half of RAW_FS in rx_stream_experiment
            # elif ((queue_idx) <  len(queue) - (len(self.REF_PREAMBLE) + 5*self.WINDOW_SIZE)): 
                possible_idx = np.array([])
                possible_doppler_shifts = np.array([])
                print("Starting packet detection.")
                # possible_idx, possible_doppler_shifts = self.pkt_detection(queue[queue_idx:], SF=self.SF, BW=self.BW, FS=self.FS, num_preamble=9)

                min_idx, max_idx = self.activity_detection(queue[queue_idx:],BW=self.BW, FS=self.FS)
                
                if (len(min_idx) > 0):  # suggests there is a packet, go thru indices
                    for SF in [7,8,9, 10, 11, 12]:
                        possible_idx_t, possible_doppler_shifts_t = self.pkt_detection(queue[queue_idx+min_idx[0]:queue_idx+max_idx[0]], SF, BW=self.BW, FS=self.FS, num_preamble=9)
                        
                        # non-zero packet detection, set the SF and proceed as usual. 
                        if len(possible_idx_t) > 0: 
                            self.SF = SF 
                            self.N = 2**SF
                            print("Current SF: ", SF, "Len: ", len(possible_idx_t))
                            # possible_idx = list(np.array(possible_idx_t) + min_idx[0])
                            possible_idx = np.array(possible_idx_t) + min_idx[0]
                            delt = np.argwhere(possible_idx < 0) 
                            possible_idx = list(possible_idx)
                            possible_doppler_shifts_t = list(possible_doppler_shifts_t)
                            for delta in delt:
                                print(delta)
                                del possible_idx[delta[0]]
                                del possible_doppler_shifts_t[delta[0]]
                            
                            print("Deleted the following invalid entries:", delt)
                            # print("Current SF: ", SF, "Len: ", possible_idx_t, len(possible_idx))
                            possible_doppler_shifts = possible_doppler_shifts_t
                # print(len(possible_idx),"ss")
                if len(possible_idx) > 0:          
                    self.PACKET_DETECTED = True
                    self.PREAMBLE_DEMOD = False
                    
                    self.possible_idx = possible_idx
                    self.possible_doppler_shifts = possible_doppler_shifts
                    if (possible_idx[0] + queue_idx > self.SNR_win):
                        self.noise_sig = queue[possible_idx[0] + queue_idx-self.SNR_win:possible_idx[0] + queue_idx]
                        # print("new noise sig:", len(self.noise_sig))
                else: 
                    self.PACKET_DETECTED = False
                    print("No packet in this queue.")
                    # queue_idx = 0 
                break;
            else:
                print("Prev queue not big enough for pkt detection.", len(queue[queue_idx:]))
                break
        
        self.LEFTOVER = queue[queue_idx:]
        self.PREV_QUEUE = []
        
    def activity_detection(self, Rx_Buffer, BW, FS): 
        '''
        Uses superDC to find energy associated with any chirp to provide an estimate location for 
        later packet detection. 
        '''
        # t_start = time.time()
        upsampling_factor = int(FS / BW)
        Pream_ind = []
        gain_history = []
        g_idx = []
        # for window in [self.superDChigh,self.superDClow]:
        for window in [self.superDClow, self.superDChigh]:
            win_size = len(window)
            loop = 0
            for off in range(3):
                offset = off * win_size // 3
                loop = Rx_Buffer.size // (win_size) - 1
                for i in range(loop):
                    temp_wind_fft = abs(
                        np.fft.fft(Rx_Buffer[(i * win_size) + offset:
                                             ((i + 1) * win_size) + offset] * np.conjugate(window), axis=0))
                    
                    noise_floor = np.mean(temp_wind_fft)
                    
                    fft_peak = np.max(temp_wind_fft)
                    
                    peak_gain = 10 * math.log10(fft_peak/noise_floor) 
                    gain_history.append(peak_gain)
                    g_idx.append(i*win_size)
                    if peak_gain > 8: 
                        Pream_ind.append(i*win_size)
      
        # # Synchronization
        Pream_ind.sort()
        # print("Activity Detection time: ", time.time() - t_start)
        if len(Pream_ind)>0:
            return [min(Pream_ind)], [max(Pream_ind) ]
        else: 
            return [], []
            
    def pkt_detection(self, Rx_Buffer, SF, BW, FS,  num_preamble):
        # t_start = time.time()
        upsampling_factor = int(FS / BW)
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
                
                b = np.argmax(temp_wind_fft)
                if len(ind_buff) >= num_preamble:
                    ind_buff = ind_buff[-(num_preamble - 1):]
                    ind_buff = np.append(ind_buff, b)                    
                else:
                    #pass
                    ind_buff = np.append(ind_buff, b)
                    #print("Overlap")
                    
                if ((sum(abs(np.diff(np.mod(ind_buff, N + 1)))) <= (num_preamble) or
                     sum(abs(np.diff(np.mod(ind_buff, N)))) <= (num_preamble) or
                     sum(abs(np.diff(np.mod(ind_buff, N - 1)))) <= (num_preamble)) and
                        ind_buff.size >= num_preamble - 1):
                    if np.sum(np.abs(Rx_Buffer[(i * upsampling_factor * N)
                                               + offset:((i + 1) * upsampling_factor * N) + offset])) != 0:
                        count = count + 1
                        Pream_ind = np.append(Pream_ind, (i - (num_preamble - 1)) * (upsampling_factor * N) + offset)
        
        # # Synchronization
        Pream_ind.sort()
        # print("Packet Detection time: ", time.time() - t_start)
        Dop_out = Pream_ind
        return list(Pream_ind), list(Dop_out)         
    def check_doppler_preamble(self, N, FS, SF, BW, num_preamble, window_size, Rx_Buffer,ind):
        '''
        Doppler shift for a packet is likely to be larger than the chirp bandwidth,
        in the especially in the narrowband case. 
        
        Dechirp successive samples and choose the FFT shift that dechirps the most energy. 
        '''
        # t_start = time.time()
        num_preamble = num_preamble - 1
        # nshifts = N*int(FS/BW) #* 2
        nshifts = N*int(FS/BW) #* 3
        upsamp = int(FS/BW)
        shifts = np.arange(-nshifts, nshifts, int(upsamp))
        shifts_small = np.arange(-int(upsamp), int(upsamp),1)
        out_preamb= 0
        possible_dop=0
        max_shift_arr = []
        max_bin_arr = []
        npt=8
        ind_old = ind
        max_val =  []
        for shift in shifts:
        # for i in [0]:
            ind_win = ((num_preamble)*window_size) + ind + shift
            win2 = np.concatenate(( np.conjugate(self.sym_2_css([0], N, SF, BW, FS)),self.sym_2_css([0], N, SF, BW, FS)))
            win1 =  Rx_Buffer[ind_win-window_size:ind_win+window_size]
            # win2 = np.conjugate(self.sym_2_css([0], N, SF, BW, FS))
            # win1 =  Rx_Buffer[ind_win-window_size:ind_win]            
            doppler_fft = np.abs(np.fft.fft(win1 * win2) )
            max_val.append(max(doppler_fft))
        
        # plt.figure(5)
        # plt.plot(shifts, max_val)
        # plt.show()
        
        st_dx = shifts[np.argmax(max_val)]
        
        max_val = []
        argmax = []
        for shift in shifts_small:
        # for i in [0]:
            ind_win = ((num_preamble)*window_size) + ind + st_dx + shift
            # win2 = np.concatenate(( np.conjugate(self.sym_2_css([0], N, SF, BW, FS)),self.sym_2_css([0], N, SF, BW, FS)))
            # win1 =  Rx_Buffer[ind_win-window_size:ind_win+window_size]
            # win2 = np.conjugate(self.sym_2_css([0], N, SF, BW, FS))
            # win1 =  Rx_Buffer[ind_win-window_size:ind_win]            
            # doppler_fft = np.abs(np.fft.fft(win1 * win2, 32*len(win2)) )
            
            
            # win1 =  Rx_Buffer[ind_win-window_size:ind_win]
            # win2 =  Rx_Buffer[ind_win:ind_win+window_size]
            # win1 =  Rx_Buffer[ind_win-2*window_size:ind_win]
            # win2 =  Rx_Buffer[ind_win:ind_win+2*window_size]
            win1 =  Rx_Buffer[ind_win-8*window_size:ind_win-4*window_size]
            win2 =  np.conjugate(Rx_Buffer[ind_win-4*window_size:ind_win])
            doppler_fft = np.abs(np.fft.fft(win1 * win2, int(len(win2))) )
            
            # win2 = np.conjugate(self.sym_2_css([0], N, SF, BW, FS))
            # win1 =  Rx_Buffer[ind_win-2*window_size:ind_win-window_size]            
            # doppler_fft2 = np.abs(np.fft.fft(win1 * win2, 32*len(win2)) )
            
            # max_val.append(max(doppler_fft) )
            # max_val.append(doppler_fft[0])
            # max_val.append(np.mean(doppler_fft[np.arange(-4,4)]))
            max_val.append(max(doppler_fft))
            argmax.append(np.argmax(doppler_fft))
        # plt.figure(5)
        # plt.plot(shifts_small, max_val)
        # plt.show()
        
        # plt.figure(6)
        # plt.plot(shifts_small,argmax)
        # plt.show()
        
        
        shift_small = shifts_small[np.argmax(max_val)]
        
        # st_dx = 0
        # shift_small = FFO_corr_UPSAMP(Rx_Buffer,st_dx+ind, SF, BW, FS, N, num_preamble-1)
        # shift_small = FFO_corr_UPSAMP(Rx_Buffer ,st_dx+ind, SF, BW, FS, N, 1)
        # shift_small = 0
        total_shifts = shift_small + st_dx #-3
        
        # plt.figure(1)
        # plt.plot(Rx_Buffer)
        # plt.axvline(ind + st_dx, color='r')
        # plt.axvline(ind + total_shifts)
        # plt.title('Blue is final shift')
        # plt.show()
        # total_shifts = shift_small + st_dx
        # if (len(max_shift_arr) > 0):
        ''' Frequency correction ... '''
        # ind = ind + shifts[np.argmax(max_shift_arr)] + 1
        # ind = ind + shifts[np.argmax(max_shift_arr)] 
        ind = ind + total_shifts
        # print("Current preamble shift: ", shifts[np.argmax(max_shift_arr)] )
        ind_win = ((num_preamble)*window_size) + ind #+ window_size
        # win2 = ( np.conjugate(self.sym_2_css([0], N, SF, BW, FS)))
        # win1 =  Rx_Buffer[ind_win-window_size:ind_win]
        win2 = ( np.conjugate(self.sym_2_css([0], N, SF, BW, FS)))
        # win1 =  Rx_Buffer[(ind_win-8*window_size):(ind_win-7*window_size)]
        win1 =  Rx_Buffer[ind:ind + window_size]
        temp_wind_fft = abs(
                np.fft.fft(win1 * win2, n=npt*len(win2), axis=0))        
        max_bin = np.argmax(temp_wind_fft)
        
        if (max_bin > ((self.N*self.UPSAMP*npt)/2)):    
            max_bin = -((self.N*self.UPSAMP*npt) - max_bin)
        freq_shift = self.FS * max_bin/(self.N*self.UPSAMP*npt)
            
        print(freq_shift, "freq shift preabmle", "Freq resolution: ", self.FS * 1/(self.N*self.UPSAMP*npt), "samp len:", len(win1))
            
        # possible_dop=int(freq_shift)
        possible_dop=float(freq_shift)
        # print("Possible frequency shift:", freq_shift)
        out_preamb=ind
        # break;
        # if possible_dop > 200: 
        # print(ind,len(Rx_Buffer),max_bin,"DOPPLER INFO")
        # plt.figure(1)
        # plt.plot(Rx_Buffer)
        # plt.axvline(ind, color='r')
        # plt.axvline(ind_old)
        # plt.figure(8)
        # plt.plot(temp_wind_fft_t[np.arange(-200,200)])
        # plt.plot(temp_wind_fft_t)
        
        plt.show()
        # print("Preamble doppler time: ", time.time() - t_start)
        return out_preamb, possible_dop         

        
    def demod_preamble(self, queue, ref_preamble, window_size, FS,freq_shift1,FD_COARSE,FD_FINE,N,SF,BW):
        '''
        Get dopppler from the preamble, and decode packet length. 
        '''   
        
        # give a few samples prior to the estimated preamble location to account for fine sampling offset     
        ind = 0 
        tt = np.arange(0, len(queue))/FS 
        num_preamble = 8
        idx_shift= FFO_corr(queue*np.exp(1j * 2 * math.pi * float(freq_shift1) * tt),ind, SF, BW, FS, N, num_preamble)
        idx_shift = 3 # TODO [REMOVE] 
        queue = queue[idx_shift:]
        print("FFO shift:", idx_shift)
        
        demod_status = False     
        
        tp = np.linspace(0, (len(ref_preamble)+2*window_size)/FS, len(ref_preamble)+2*window_size)
        
        
        freq_shift = freq_shift1    
        ts = np.linspace(0, window_size/FS, window_size)
        
        # check to make sure we can demodulate the preamble. If not, something went wrong with our sync and we need to try another pt 
        
        doppler_sample = np.ones(window_size)
        for symdx in range(0,8):
            # print(self.symbol_demod_sig(self.PREV_QUEUE[self.WINDOW_SIZE*(symdx):self.WINDOW_SIZE*(symdx+1)]))
            # doppler = (np.arange(0, window_size) / FS ) * 5
            # sym_win = queue[window_size*(symdx):window_size*(symdx+1)] * np.exp(1j * 2 * math.pi * doppler * ts)
            
            sym_win = queue[window_size*(symdx):window_size*(symdx+1)] * np.exp(1j * 2 * math.pi * freq_shift * ts)
            # sym_win = queue[window_size*(symdx):window_size*(symdx+1)] * np.exp(1j * 2 * math.pi * freq_shift * ts) * np.conjugate(doppler_sample)
            dem_sym = self.symbol_demod_sig(sym_win) 
            freq_shift2 = FFO_corr_sym(sym_win, self.SF, self.BW, self.FS, self.N, 0)
            
            freq_shift = freq_shift + freq_shift2
            print("New f shift: " , freq_shift2, "Total f shift: ",freq_shift)
            
            if (dem_sym != 0):
                print("Sync or doppler correction is incorrect.",self.symbol_demod_sig(sym_win))
                print("Sym dx: ", symdx)
                demod_status = False
                return demod_status, 0,0,0  
                        
        for downchirp in range(0,2): 
            # detect the downchirps (sync word) 
            # sym_win = queue[window_size*(downchirp+8):window_size*(downchirp+9)] * np.exp(1j * 2 * math.pi * freq_shift * ts)
            # sym_win = queue[window_size*8 + (window_size*downchirp):window_size*(downchirp+8)+(window_size*2*downchirp)] * np.exp(1j * 2 * math.pi * freq_shift * ts)
            st_dx = window_size*8 + downchirp*window_size 
            end_dx = st_dx+window_size
            sym_win = queue[st_dx:end_dx] * np.exp(1j * 2 * math.pi * freq_shift * ts)
            
            # plt.figure(1)
            # plt.specgram(sym_win)
            # plt.show()
            
            # dem_sym = self.symbol_demod_sig_dc(sym_win)
            dem_sym = self.symbol_demod_sig(np.conjugate(sym_win))            
            # print("Preamble symbol: ", dem_sym)
            # freq_shift2 = FFO_corr_sym_conj(sym_win, self.SF, self.BW, self.FS, self.N, 0)
            # freq_shift2 = 0
            freq_shift2 = -FFO_corr_sym(np.conjugate(sym_win), self.SF, self.BW, self.FS, self.N, 0)
            freq_shift = freq_shift + freq_shift2
            print("New f shift [coj]: " , freq_shift2, dem_sym, "Total f shift: ",freq_shift)
            
            if (dem_sym != 0):
                print("Sync or doppler correction is incorrect. [DC]",self.symbol_demod_sig(sym_win))
                print("Sym dx: ", downchirp, dem_sym)
                demod_status = False
                return demod_status, 0,0,0  
            else:
                print("Passed.")        
        t = ts = np.linspace(0, len(queue)/FS, len(queue))
        # queue = queue * np.exp(1j * 2 * math.pi * freq_shift * t)   
        self.OUTPUT_PKT = []
        end_preamb_dx = 2*window_size+len(ref_preamble)
        
        # freq_shift = freq_shift + freq_shift2        
        pkt_len_2 = self.symbol_demod_sig(queue[end_preamb_dx:end_preamb_dx+window_size]* np.exp(1j * 2 * math.pi * freq_shift * t[:window_size]) )  
            
        freq_shift2 = FFO_corr_sym(queue[end_preamb_dx:end_preamb_dx+window_size]* np.exp(1j * 2 * math.pi * freq_shift * t[:window_size]), self.SF, self.BW, self.FS, self.N, pkt_len_2)
        # print("-New f shift: " , freq_shift2)
        freq_shift = freq_shift+ freq_shift2
            
        pkt_len_1 = self.symbol_demod_sig(queue[end_preamb_dx+window_size:end_preamb_dx+2*window_size]* np.exp(1j * 2 * math.pi * freq_shift * t[:window_size]))
        
        freq_shift2 = FFO_corr_sym(queue[end_preamb_dx+window_size:end_preamb_dx+2*window_size]* np.exp(1j * 2 * math.pi * freq_shift * t[:window_size]), self.SF, self.BW, self.FS, self.N, pkt_len_1)
        # print("-New f shift: " , freq_shift2)
        freq_shift = freq_shift+ freq_shift2       
                        
        #[TODO] : What if other candidate IDXs lie within this, and this sync is just a bit off? Need to not remove the pts then
        packet_len =  pkt_len_2 + N*pkt_len_1
        if (packet_len > 1000):
            print(" ===== found pkt length too long ", packet_len, pkt_len_2, pkt_len_1)
            pass
        else: 
            demod_status = True    
        # print("Packet length:", packet_len)    
        return demod_status, packet_len, freq_shift, idx_shift        
    # def demod_preamble(self, queue, ref_preamble, window_size, FS,freq_shift1,FD_COARSE,FD_FINE,N,SF,BW):
        # '''
        # Get dopppler from the preamble, and decode packet length. 
        # '''   
        # # plt.figure(1)
        # # plt.plot(queue)
        # # plt.show()
        # print("Demodulating the preamble.")
        # # give a few samples prior to the estimated preamble location to account for fine sampling offset     
        # ind = 0 
        # tt = np.arange(0, len(queue))/FS 
        # num_preamble = 8
        # idx_shift= FFO_corr(queue*np.exp(1j * 2 * math.pi * float(freq_shift1) * tt),ind, SF, BW, FS, N, num_preamble)
        # queue = queue[idx_shift:]
        # print("FFO shift:", idx_shift)
        
        # demod_status = False     
        
        # tp = np.linspace(0, (len(ref_preamble)+2*window_size)/FS, len(ref_preamble)+2*window_size)
        
        # # preamb_with_sync = np.concatenate((ref_preamble, np.conjugate(self.sym_2_css([0,0], N, SF,BW,FS))))
        # # preamb_candidate = queue[:len(ref_preamble)+window_size*2]
        # # preamb_candidate = preamb_candidate * np.exp(1j * 2 * math.pi * freq_shift1 * tp)
        
        # freq_shift = freq_shift1    
        # ts = np.linspace(0, window_size/FS, window_size)
        # # [TODO] why was tihs here?
        # # freq_shift2 = FFO_corr_sym(queue[window_size*(0):window_size*(0+1)] * np.exp(1j * 2 * math.pi * freq_shift * ts), self.SF, self.BW, self.FS, self.N, 0)
        # # # print("Init fine f shift: " , freq_shift2)
        # # freq_shift = freq_shift + freq_shift2
        
        # # check to make sure we can demodulate the preamble. If not, something went wrong with our sync and we need to try another pt 
        
        # for symdx in range(0,8):
            # # print(self.symbol_demod_sig(self.PREV_QUEUE[self.WINDOW_SIZE*(symdx):self.WINDOW_SIZE*(symdx+1)]))
            # # sym_win = queue[window_size*(symdx):window_size*(symdx+1)] * np.exp(1j * 2 * math.pi * freq_shift1 * ts)
            # # sym_win = queue[window_size*(symdx):window_size*(symdx+1)] * np.exp(1j * 2 * math.pi * freq_shift * ts)
            # sym_win = queue[window_size*(symdx):window_size*(symdx+1)] * np.exp(1j * 2 * math.pi * freq_shift * ts)
            
            # dem_sym = self.symbol_demod_sig(sym_win) 
            # # print("Preamble symbol: ", dem_sym)
            # freq_shift2 = FFO_corr_sym(sym_win, self.SF, self.BW, self.FS, self.N, 0)
            # # freq_shift2 = 0
            # print("New f shift: " , freq_shift2)
            # freq_shift = freq_shift + freq_shift2
            # if (dem_sym != 0):
                # print("Sync or doppler correction is incorrect.",self.symbol_demod_sig(sym_win))
                # print("Sym dx: ", symdx)
                # demod_status = False
                # return demod_status, 0,0,0             
        
        # for downchirp in range(0,2): 
            # # detect the downchirps (sync word) 
            # sym_win = queue[window_size*(downchirp+8):window_size*(downchirp+9)] * np.exp(1j * 2 * math.pi * freq_shift * ts)
            
            # # plt.figure(1)
            # # plt.specgram(sym_win)
            # # plt.show()
            
            # # dem_sym = self.symbol_demod_sig_dc(sym_win)
            # dem_sym = self.symbol_demod_sig(np.conjugate(sym_win))            
            # # print("Preamble symbol: ", dem_sym)
            # # freq_shift2 = FFO_corr_sym_conj(sym_win, self.SF, self.BW, self.FS, self.N, 0)
            # # freq_shift2 = 0
            # freq_shift2 = FFO_corr_sym(np.conjugate(sym_win), self.SF, self.BW, self.FS, self.N, 0)
            
            # # print("New f shift [coj]: " , freq_shift2)
            # # freq_shift = freq_shift + freq_shift2
            # if (dem_sym != 0):
                # print("Sync or doppler correction is incorrect. [DC]",self.symbol_demod_sig(sym_win))
                # print("Sym dx: ", downchirp, dem_sym)
                # demod_status = False
                # return demod_status, 0,0,0  
            # else:
                # print("Passed.")
        # t = ts = np.linspace(0, len(queue)/FS, len(queue))
        # # queue = queue * np.exp(1j * 2 * math.pi * freq_shift * t)   
        # self.OUTPUT_PKT = []
        # end_preamb_dx = 2*window_size+len(ref_preamble)
        
        # # doppler hidden by the two upchirps here ... 
        # # the down chirps also accumulate frequency shifts 
        # # freq_shift2 = FFO_corr_sym_conj(queue[len(ref_preamble):(len(ref_preamble))+window_size]* np.exp(1j * 2 * math.pi * freq_shift * t[:window_size]), self.SF, self.BW, self.FS, self.N, 0)
        # # print("-!New f shift: " , freq_shift2)
        # # freq_shift = freq_shift+ freq_shift2
        
        # # freq_shift2 = FFO_corr_sym_conj(queue[(len(ref_preamble))+window_size:(len(ref_preamble))+2*window_size]* np.exp(1j * 2 * math.pi * freq_shift * t[:window_size]), self.SF, self.BW, self.FS, self.N, 0)
        # # print("-!New f shift: " , freq_shift2)
        # # freq_shift = freq_shift+ freq_shift2

        
        # # freq_shift = freq_shift + freq_shift2        
        # pkt_len_2 = self.symbol_demod_sig(queue[end_preamb_dx:end_preamb_dx+window_size]* np.exp(1j * 2 * math.pi * freq_shift * t[:window_size]) )  
            
        # freq_shift2 = FFO_corr_sym(queue[end_preamb_dx:end_preamb_dx+window_size]* np.exp(1j * 2 * math.pi * freq_shift * t[:window_size]), self.SF, self.BW, self.FS, self.N, pkt_len_2)
        # print("-New f shift: " , freq_shift2)
        # freq_shift = freq_shift+ freq_shift2
            
        # pkt_len_1 = self.symbol_demod_sig(queue[end_preamb_dx+window_size:end_preamb_dx+2*window_size]* np.exp(1j * 2 * math.pi * freq_shift * t[:window_size]))
        
        # freq_shift2 = FFO_corr_sym(queue[end_preamb_dx+window_size:end_preamb_dx+2*window_size]* np.exp(1j * 2 * math.pi * freq_shift * t[:window_size]), self.SF, self.BW, self.FS, self.N, pkt_len_1)
        # print("-New f shift: " , freq_shift2)
        # freq_shift = freq_shift+ freq_shift2       
                        
        # print("Final preamble frequency shift: ", freq_shift)
        # #[TODO] : What if other candidate IDXs lie within this, and this sync is just a bit off? Need to not remove the pts then
        # packet_len =  pkt_len_2 + N*pkt_len_1
        # if (packet_len > 1000):
            # print(" ===== found pkt length too long ", packet_len)
            # print("pkt_len_2: ",pkt_len_2)
            # print("pkt_len_1: ",pkt_len_1)
            # pass
        # else: 
            # demod_status = True    
        # # print("Packet length:", packet_len)    
        # return demod_status, packet_len, freq_shift, idx_shift
    
    def calc_SNR(self, signal_sig, noise_sig):
        '''
        Calculates the SNR based on a measured signal with/without signal of interest.
        
        Should be downsampled to the appropriate BW before computation. 
        '''    

        noise_power = np.sum(np.abs(noise_sig)**2) / (len(noise_sig))
        
        signal_power = np.sum(np.abs(signal_sig)**2) / len(signal_sig)           
             
        current_SNR = 20 * np.log10( (signal_power-noise_power) / noise_power) 
        return current_SNR
        
    def symbol_demod_sig(self,sig_tmp):
        '''
        Demodulates a CSS symbol and returns the frequency bin 
        at which the symbol appears.
        '''        
        
        trans_upchirp = np.conjugate(self.UPCHIRP)
        
        sig_tmp = scipy.signal.decimate(sig_tmp, self.UPSAMP)
        trans_upchirp = scipy.signal.decimate(trans_upchirp, self.UPSAMP)
        
        # trans_upchirp = trans_upchirp[::10]
        dechirped = sig_tmp * trans_upchirp
                
        dechirped = np.squeeze(dechirped)

        # dechirped = dechirped[::10]

        # data_fft = abs(np.fft.fft(dechirped)).transpose()
        data_fft = abs(np.fft.fft(dechirped))
        # dechirped = np.concatenate((data_fft[:int(self.N/2)], \
             # data_fft[-int(self.N/2):]))
        dechirped = data_fft
    
        freq_bin = np.argmax(dechirped)    
                                
        return freq_bin
    def symbol_demod_sig_dc(self,sig_tmp):
        '''
        Demodulates a CSS symbol and returns the frequency bin 
        at which the symbol appears.
        '''        
        
        trans_upchirp = np.conjugate(self.UPCHIRP)
        # trans_upchirp = self.UPCHIRP
        
        sig_tmp = np.conjugate(sig_tmp)
        sig_tmp = scipy.signal.decimate(sig_tmp, self.UPSAMP)
        trans_upchirp = scipy.signal.decimate(trans_upchirp, self.UPSAMP)
        
        # trans_upchirp = trans_upchirp[::10]
        dechirped = sig_tmp * trans_upchirp
        
            
        dechirped = np.squeeze(dechirped)

        # dechirped = dechirped[::10]

        # data_fft = abs(np.fft.fft(dechirped)).transpose()
        data_fft = abs(np.fft.fft(dechirped))
        
        plt.figure(67)
        plt.plot(data_fft)
        plt.show()
        # dechirped = np.concatenate((data_fft[:int(self.N/2)], \
             # data_fft[-int(self.N/2):]))
        dechirped = data_fft
    
        freq_bin = np.argmax(dechirped)    
                                
        return freq_bin
    def create_upchirp(self):
        '''
        Create an upchirp which is used to demodulate the symbols 
        '''
        return self.sym_2_css([0], self.N, self.SF, self.BW, self.FS)
    
    def genSuperDC(self, SF, BW, FS): 
        '''
        Generate a super downchirp consisting of SF, SF-1, and SF-2 chirps 
        '''
        # print(len(np.tile(self.sym_2_css([0], N, SF, BW, FS),1)),len(np.tile(self.sym_2_css([0], N, SF-1, BW, FS),2)))
        return np.tile(self.sym_2_css([0], 2**(SF), SF, BW, FS),1) + \
        np.tile(self.sym_2_css([0], 2**(SF-1), SF-1, BW, FS),2) + \
        np.tile(self.sym_2_css([0], 2**(SF-2), SF-2, BW, FS),4)
    
    def sym_2_css(self, symbol, N, SF, BW, Fs): 
        '''
        sym_2_css :: Modulates a symbol between (0,N-1) to a chirp  
        
        N = 2^SF == does not represent upsampling factor (?)
        SF 
        BW = Bandwidth sweep of the chirp -BW/2 to BW/2 
        Fs = Sampling rate
        '''
        
        if (N != 2**SF): 
            return
        
        N = 2**SF
        sym_out = []
        
        # Fs/Bw is upsamp...
        spsym = int(Fs/BW*N) # symbols defined by their freq offset at the start 
        # spsym = N*10
        
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
        if (len(self.GND_TRUTH_PKT) > 1): 
            # print("!!!!!Received packet at time: ", self.CURR_POSIX_TIME)
            with open(self.OUTFILE+'_'+str(self.TOTAL_PKT_CNT) + '.pkl', 'wb') as f: 
                pickle.dump([self.GND_TRUTH_PKT,  self.OUTPUT_PKT, self.TOTAL_PKT_CNT,self.PKT_SNR, self.CURR_POSIX_TIME],f)
            print("Ground Truth:",self.GND_TRUTH_PKT)
            print("Output:      ", self.OUTPUT_PKT)
            if (len(self.OUTPUT_PKT) == len(self.GND_TRUTH_PKT)):
                nwrong = np.count_nonzero(np.abs(np.subtract(self.OUTPUT_PKT , self.GND_TRUTH_PKT)))
                print("num wrong: ",nwrong)    
                if (nwrong > 0):
                    print('\007')                
            else: 
                print("Packet length does not agree." )
                print('\007')  
            print("TOTAL PACKET COUNT: ", self.TOTAL_PKT_CNT)
        else: 
            print("Output:      ", self.OUTPUT_PKT)
            with open(self.OUTFILE+'_'+str(self.TOTAL_PKT_CNT) + '.pkl', 'wb') as f: 
                pickle.dump([self.GND_TRUTH_PKT,  self.OUTPUT_PKT, self.TOTAL_PKT_CNT,self.PKT_SNR, self.CURR_POSIX_TIME],f)
            
    # def get_doppler_preamble(self, reference, rx_data, doppler_freq,step,t): 
        # '''
        # Corrects doppler shift based on the known preamble 
        # '''       
        # f_d = np.arange(-doppler_freq, doppler_freq, step)
        # xcorr_arr = []
        # for f in f_d: 
            # t = np.linspace(0, len(rx_data)/self.FS, len(rx_data))
            # data_shifted = rx_data * np.exp(1j * 2 * math.pi * f * t)   
            # xcorr_val = np.squeeze(abs(np.correlate(data_shifted, reference)))            
            # xcorr_arr.append(xcorr_val)
                  
        # freq_bin_shift = f_d[np.argmax(xcorr_arr)]
        
        # return freq_bin_shift
        
    # def get_doppler_symbol(self, reference, rx_data, doppler_freq,step,t): 
        # '''
        # Corrects doppler shift based on the known preamble 
        # '''       
        # f_d = np.arange(-doppler_freq, doppler_freq, step)
        # xcorr_arr = []
        # for f in f_d: 
            # t = np.linspace(0, len(rx_data)/self.FS, len(rx_data))
            # data_shifted = rx_data * np.exp(1j * 2 * math.pi * f * t)   
            # xcorr_val = np.squeeze(abs(np.correlate(data_shifted, reference)))            
            # xcorr_arr.append(xcorr_val)
                  
        # freq_bin_shift = f_d[np.argmax(xcorr_arr)]
        
        # return freq_bin_shift