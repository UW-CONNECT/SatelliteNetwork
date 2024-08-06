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
import statistics as st # from APD
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
        self.prev_sym = 0
        self.curr_sym = 0
                
        self.doppler_slope = 0
        
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
                and (max(self.possible_idx)) < (len(queue)-(len(self.REF_PREAMBLE) + self.WINDOW_SIZE*10))):
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
                    # print("Out of doppler preamble.", len(self.possible_idx), "PDS")
                    # end_preamble_dx = possible_idx + len(self.REF_PREAMBLE) + self.WINDOW_SIZE*4 #todo [sub in automatic length]
                    end_preamble_dx = possible_idx + self.WINDOW_SIZE*12 #todo [sub in automatic length]
                    
                    if (end_preamble_dx > len(queue) or possible_idx < 0): # 6/5/2024 change this possible 
                        print("Breaking")
                        possible_idx = queue_idx # may need to re-do pkt detection?
                        break;
                    
                    # plt.figure(18)
                    # plt.plot(queue)
                    # plt.axvline(possible_idx)
                    # plt.title("Sync check")
                    # plt.show()
                    
                    ''' FAKE DOPPLER SHIFT ADDITION '''  
                    
                    # plt.figure(1)
                    # plt.plot(fake_dop)
                    # plt.title('Fake doppler')
                    # plt.show()
                    # queue_tmp = queue[possible_idx:end_preamble_dx+14] 
                    queue_tmp = queue 
                    tt = np.arange(0, len(queue_tmp[possible_idx:])) / self.FS 
                    doppler_slope = 50
                    fake_dop = tt * doppler_slope
                    # fake_dop_vec = np.exp(1j * 2 * math.pi * doppler_slope * tt**2/2)
                    fake_dop_vec = 1
                    # print("Expected doppler value: ", doppler_slope * tt[-1])
                    queue_tmp[possible_idx:] = queue_tmp[possible_idx:]* fake_dop_vec
                    # Freq = (np.angle(fake_dop_vec[1:] * np.conjugate(fake_dop_vec[0:-1])));
                    # Freq = Freq / (2*math.pi) * (self.FS);      
                    
                    # plt.figure(1)
                    # plt.plot(Freq)
                    # plt.show()

                    # print(possible_idx, end_preamble_dx, len(queue_tmp))
                    self.PREAMBLE_DEMOD, self.PACKET_LEN, freq_shift,idx_shift = self.demod_preamble(queue_tmp[possible_idx:end_preamble_dx+14], self.REF_PREAMBLE, self.WINDOW_SIZE, self.FS, -possible_dop, self.FD_COARSE, self.FD_FINE,self.N,self.SF,self.BW)
                    # print("Attempting preamble demod. Possible idx:", possible_idx, "possible doppler:", possible_dop)
                    # self.PREAMBLE_DEMOD, self.PACKET_LEN, freq_shift,idx_shift = self.demod_preamble(queue[possible_idx-3:end_preamble_dx+14], self.REF_PREAMBLE, self.WINDOW_SIZE, self.FS, -possible_dop, self.FD_COARSE, self.FD_FINE,self.N,self.SF,self.BW)
                    # self.PREAMBLE_DEMOD, self.PACKET_LEN, freq_shift,idx_shift = self.demod_preamble(queue[possible_idx:end_preamble_dx], self.REF_PREAMBLE, self.WINDOW_SIZE, self.FS, -possible_dop, self.FD_COARSE, self.FD_FINE,self.N,self.SF,self.BW)          # [TODO] 7/1/2024 remove
                    # demod_status, packet_len, freq_shift
                    if (self.PREAMBLE_DEMOD):
                        # [TODO] Remove 
                        queue = queue_tmp
                        
                        # '''
                        # Apply doppler slope to the whole packet, so that all following doppler corrections
                        # may be regarded as a constant doppler shift. 
                        # '''
                        # obvs_idx = possible_idx + len(self.REF_PREAMBLE) + self.WINDOW_SIZE*4 +idx_shift-3
                        # tt = np.arange(0, len(queue[obvs_idx:])) / self.FS 
                        # # doppler_slope = 50
                        # fake_dop = tt * (self.doppler_slope/tt[self.WINDOW_SIZE])
                        # print("DOPPLER SLOPE HZ/s ", (self.doppler_slope/tt[self.WINDOW_SIZE]))
                        # doppler_slope_component =  np.exp(1j * 2 * math.pi * doppler_slope * tt**2/2)
                        # doppler_offset = np.exp(1j * 2 * math.pi * freq_shift * tt) 
                        # queue[obvs_idx:] = queue[obvs_idx:]*doppler_slope_component*doppler_offset
                                            
                        self.END_COUNTER = 0 
                        self.PACKETS_DECODED = 0 
                        # self.doppler_cum = self.doppler_cum + freq_shift
                        self.doppler_cum = freq_shift
                        # self.doppler_cum = 0
                        
                        # print("Doppler after preamble demod: ", self.doppler_cum)
                        # calculate SNR if packet detected 
                        # ds_factor = int(self.WINDOW_SIZE / (self.N * 10))  
                        ds_factor = self.UPSAMP
                        # print("Noise sig len:", len(self.noise_sig))
                        # noise_sig = scipy.signal.decimate(self.noise_sig, ds_factor)
                        noise_sig = self.noise_sig[::ds_factor]
                        
                        # print("Noise sig len:", len(noise_sig))
                        # noise_sig = self.noise_sig
                        SNR_win = self.SNR_win
                        # signal_sig = queue[possible_idx:possible_idx+SNR_win]
                        tt = np.arange(0, SNR_win)/self.FS
                        signal_sig = queue[possible_idx:possible_idx+SNR_win]* np.exp(1j * 2 * math.pi * self.doppler_cum * tt).transpose()
                        
                        # print("Peak gain SNR: ", 2*self.peak_gain_calculation(np.abs(np.fft.fft(signal_sig[:self.WINDOW_SIZE]*np.conjugate(self.UPCHIRP  )))))
                        
                        # plt.figure(1)
                        # plt.plot(np.abs(np.fft.fft(signal_sig[:self.WINDOW_SIZE]*np.conjugate(self.UPCHIRP  ))))
                        # plt.show()
                        
                        # signal_sig = scipy.signal.decimate(signal_sig, ds_factor)
                        signal_sig = signal_sig[::ds_factor]                        

                        self.PKT_SNR = self.calc_SNR(signal_sig, noise_sig)
                        # print("SNR For this packet:",self.PKT_SNR, "Init freq shift: ", freq_shift, "Packet len:", self.PACKET_LEN)
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
                    # tt = np.arange(0, 2*self.WINDOW_SIZE)/self.FS
                    tt = np.arange(0, self.WINDOW_SIZE)/self.FS
                    
                    sym = self.symbol_demod_sig(queue[queue_idx:queue_idx+self.WINDOW_SIZE]* np.exp(1j * 2 * math.pi * self.doppler_cum * tt[:self.WINDOW_SIZE]))

                    # get a new fine-grained doppler estimate  sym_2_css(self, symbol, N, SF, BW, Fs)
                    reference = self.sym_2_css([sym], self.N, self.SF, self.BW, self.FS)
                    t_queue = np.arange(0, len(queue))/self.FS
                    freq_lin = self.BW / (self.N)
                    
                    
                    f_shift_new = FFO_corr_sym(queue[queue_idx:queue_idx+self.WINDOW_SIZE]* np.exp(1j * 2 * math.pi * self.doppler_cum * tt), self.SF, self.BW, self.FS, self.N, sym)
                    # self.curr_sym = sym 
                    # f_shift_new = doppler_slope_offset(queue[queue_idx-self.WINDOW_SIZE:queue_idx+self.WINDOW_SIZE]* np.exp(1j * 2 * math.pi * self.doppler_cum * tt), self.SF, self.BW, self.FS, self.N, self.prev_sym, self.curr_sym)     #[TESTING]               # f_shift_new=0 #[TODO] Remove.
                    # self.prev_sym = self.curr_sym 
                    # self.doppler_cum = self.doppler_cum #+ f_shift_new
                    # print("New symbol shift: ", f_shift_new)
                    self.doppler_cum = self.doppler_cum + self.doppler_slope + f_shift_new
                    # print("New symbol shift: ", f_shift_new, "Total dop: ", self.doppler_cum)
                    # self.doppler_cum = self.doppler_cum + f_shift_new
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
                        print("Doppler slope: ",self.doppler_slope)      
                    self.END_COUNTER = self.END_COUNTER + 1 
                    #self.PREAMBLE_DEMOD = False
                else: 
                    print("Ending DOP: ",self.doppler_cum)  
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
                # print("Starting packet detection.")
                # possible_idx, possible_doppler_shifts = self.pkt_detection(queue[queue_idx:], SF=self.SF, BW=self.BW, FS=self.FS, num_preamble=9)
                print("Starting packet detection.")
                # possible_idx, possible_doppler_shifts = self.pkt_detection(queue[queue_idx:], SF=self.SF, BW=self.BW, FS=self.FS, num_preamble=9)

                min_idx, max_idx = self.activity_detection(queue[queue_idx:],BW=self.BW, FS=self.FS)
                
                if (len(min_idx) > 0):  # suggests there is a packet, go thru indices
                    possible_SF=[7,8,9, 10, 11, 12]
                    past_max = 0
                    for SF in possible_SF:
                        possible_idx_t, possible_doppler_shifts_t = self.pkt_detection(queue[queue_idx+min_idx[0]:queue_idx+max_idx[0]], SF, BW=self.BW, FS=self.FS, num_preamble=9)
                                            
                    
                        # non-zero packet detection, set the SF and proceed as usual. 
                        if len(possible_idx_t) > past_max: 
                            self.SF = SF 
                            self.N = 2**SF
                            self.WINDOW_SIZE = self.N * self.UPSAMP
                            self.UPCHIRP = self.create_upchirp()
                            self.REF_PREAMBLE = self.sym_2_css([0,0,0,0,0,0,0,0], self.N, self.SF, self.BW, self.FS) # [TODO] most locations where this is used can just use 8*win_Size
                            print("Current SF: ", SF, "Len: ", len(possible_idx_t))
                            # possible_idx = list(np.array(possible_idx_t) + min_idx[0])
                            possible_idx = np.array(possible_idx_t) + min_idx[0]
                            delt = np.argwhere(possible_idx < 0) 
                            possible_idx = list(possible_idx)
                            possible_doppler_shifts_t = list(possible_doppler_shifts_t)
                            ctr = 0 
                            for delta in delt:
                                print(delta)
                                del possible_idx[delta[0] - ctr ]
                                del possible_doppler_shifts_t[delta[0] - ctr]
                                ctr = ctr + 1 
                            print("Deleted the following invalid entries:", delt)
                            # print("Current SF: ", SF, "Len: ", possible_idx_t, len(possible_idx))
                            possible_doppler_shifts = possible_doppler_shifts_t
                        # else: 
                            past_max = len(possible_idx_t)
                # # nfl = 0 # noise floor ? 
                # # uplink_wind = self.Active_Sess_SFAgn(queue[queue_idx:], int(self.FS/self.BW), 7,2,nfl)
                # # uplink_values = uplink_wind[-1][3:len(uplink_wind[-1])]
                # # print("SFAGN thing: ", uplink_wind[-1][2:len(uplink_wind[-1])])
                # # SFcnt = 7+ np.argmax(uplink_values)
                # # print("Predicted SF: ", SFcnt)
                # # # print("Uplink wind VS Rx Buffer size: ", len(uplink_wind[0]), len(Rx_Buffer))
                # # print(uplink_wind)
                # # print(f'\nFound {len(uplink_wind)} active Periods')

                # # # min_idx, max_idx = self.activity_detection(queue[queue_idx:],BW=self.BW, FS=self.FS)
                
                # # # if (len(min_idx) > 0):  # suggests there is a packet, go thru indices
                    # # # possible_SF=[7,8,9, 10, 11, 12]
                    # # # past_max = 0
                    # # # for SF in possible_SF:
                # # # possible_idx_t, possible_doppler_shifts_t = self.pkt_detection(queue[queue_idx+min_idx[0]:queue_idx+max_idx[0]], SF, BW=self.BW, FS=self.FS, num_preamble=9)
                # # possible_idx_t, possible_doppler_shifts_t = self.pkt_detection(queue[queue_idx:], SFcnt, BW=self.BW, FS=self.FS, num_preamble=9)
                                           
                    
                # # # non-zero packet detection, set the SF and proceed as usual. 
                # # # if len(possible_idx_t) > past_max: 
                # # self.SF = SFcnt 
                # # self.N = 2**self.SF
                # # self.WINDOW_SIZE = self.N * self.UPSAMP
                # # self.UPCHIRP = self.create_upchirp()
                # # print("Current SF: ", self.SF, "Len: ", len(possible_idx_t))
                # # # possible_idx = list(np.array(possible_idx_t) + min_idx[0])
                # # possible_idx = np.array(possible_idx_t) #+ min_idx[0]
                # # delt = np.argwhere(possible_idx < 0) 
                # # possible_idx = list(possible_idx)
                # # possible_doppler_shifts_t = list(possible_doppler_shifts_t)
                # # ctr = 0 
                # # for delta in delt:
                    # # print(delta)
                    # # del possible_idx[delta[0] - ctr ]
                    # # del possible_doppler_shifts_t[delta[0] - ctr]
                    # # ctr = ctr + 1 
                # # # print("Deleted the following invalid entries:", delt)
                # # # print("Current SF: ", SF, "Len: ", possible_idx_t, len(possible_idx))
                # # possible_doppler_shifts = possible_doppler_shifts_t
                # else: 
                    # past_max = len(possible_idx_t)
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
                break;  #[TODO] Remove this [7/24/2024]
            else:
                print("Prev queue not big enough for pkt detection.", len(queue[queue_idx:]))
                break
        
        self.LEFTOVER = queue[queue_idx:]
        self.PREV_QUEUE = []
    def Active_Sess_SFAgn(self, x_1, UPSAMPLE_FACTOR, minSF,win_jump_factor,nfl):
        '''
        See BYOG/Active_Period_Detector.py for appropriate values. 
        win_jump_factor = 2
        '''
        upsampling_factor = UPSAMPLE_FACTOR
        windsz = 2 ** minSF
        win_jump = int(np.floor(windsz * upsampling_factor / win_jump_factor))

        peak_gain = []
        uplink_wind = []
        n = []
        p = []
        last_wind = 0
        num_accum = 32
        # front_buf = 58 * self.win_jump_factor
        # back_buf = 78 * self.win_jump_factor
        front_buf = 100 * win_jump_factor
        back_buf = 100 * win_jump_factor
        mov_thresh_wind = int(np.floor((1000 * win_jump_factor) / 32))
        mov_thresh_rec = []
        # print(f"nfl Val = {self.nfl}")
        upLim = nfl  # db  [JMS] I am guessing this is noise floor. How to get this estimate? Calibration step?
        # print(upLim)
        # upLim = 0.4  # db
        c = 1
        t_n = 0
        t_p = 0
        x = []
        lst_ind = 0
        ctr = 0
        dis = 0
        ind_arr = []
        big_ind_arr = []
        t_u = 4.5  # 1.06
        t_v = 4.5  # 1.2
        idx = np.concatenate([np.array(range(0, int(windsz / 2))), np.array(
            range(int(windsz / 2 + (upsampling_factor - 1) * windsz), int((upsampling_factor) * windsz)))])

        # [JMS] Not sure what this is doing. is lst_chunk to be added from the buffeR?
        # # if len(self.lst_chunk) != 0:
            # # x_1 = np.append(self.lst_chunk[len(self.lst_chunk) - 16384: len(self.lst_chunk)],
                            # # x_1)

        x_1_len = len(x_1)

        # for i in range(1, int(np.floor(x_1_len / win_jump) - (win_jump_factor * (win_jump_factor - 1)))):
        # print(f"Heyyyyyyyyyyyyyyyyyy!! : {x_1_len / win_jump}")
        # for i in range(1, int(np.floor(x_1_len / win_jump))):
        for i in range(1, int(np.floor(x_1_len / win_jump) - (win_jump_factor * (win_jump_factor - 1)))):
            wind_fft = np.abs(np.fft.fft(np.multiply(x_1[(i - 1) * win_jump: (i - 1) * win_jump + (windsz * upsampling_factor)], np.conj(x_1[(i - 1) * win_jump + (windsz * upsampling_factor): (i - 1) * win_jump + (2 * windsz * upsampling_factor)]))))
            w_f = wind_fft[idx]
            ind = np.argmax(w_f)
            id = np.setdiff1d(range(0, len(w_f)), ind)
            n_f = np.mean(w_f[id])

            noise_floor = np.mean(wind_fft)
            fft_peak = max(wind_fft)

            # if len(peak_gain) < 200:
            #     upLim = t_u * np.mean(peak_gain)

            if np.mod(c, num_accum) == 0:
                n.append(t_n)
                p.append(t_p)
                peak_gain.append(t_p)
                x.append(i)

######################################################################################################################################
                # if ((i + 1) / num_accum) > mov_thresh_wind:
                #     mov_thresh = t_v * np.mean(peak_gain[len(peak_gain) - mov_thresh_wind:len(peak_gain)])
                #     if mov_thresh > upLim:
                #         mov_thresh = upLim
                # else:
                #     mov_thresh = t_v * np.mean(peak_gain)
                #     if mov_thresh > upLim:
                #         mov_thresh = upLim
######################################################################################################################################
                mov_thresh = upLim
######################################################################################################################################
                mov_thresh_rec.append(mov_thresh)

                ## Add alternative to energy thresholding
                # if len(ind_arr) > 0:
                #     D_1 = 0
                #     D_65 = 0
                #     D_97 = 0
                #     D_113 = 0
                #     D_121 = 0
                #     D_125 = 0
                #     # temp = np.array(big_ind_arr).reshape(1, -1)
                #     for o in ind_arr:
                #         if o == 0:
                #             D_1 += 1
                #         elif o == 64:
                #             D_65 += 1
                #         elif o == 96:
                #             D_97 += 1
                #         elif o == 112:
                #             D_113 += 1
                #         elif o == 120:
                #             D_121 += 1
                #         elif o == 124:
                #             D_125 += 1
                # num = 10
                # if D_1 >= num or D_65 >= num or D_97 >= num or D_113 >= num or D_121 >= num or D_125 >= num:
                #########################################

                if peak_gain[-1] >= mov_thresh:
                    big_ind_arr.append(ind_arr)
                    if i - num_accum >= last_wind:
                        if i - back_buf < 1:
                            # print('Touched this Corner Case\n')
                            uplink_wind.append([1, i + front_buf, 0, 0, 0, 0, 0, 0, 0])
                            # uplink_wind.append([i - back_buf, i + front_buf, 0, 0, 0, 0, 0, 0, 0])
                            big_ind_arr = []
                            big_ind_arr.append(ind_arr)
                        else:
                            uplink_wind.append([i - back_buf, i + front_buf, 0, 0, 0, 0, 0, 0, 0])
                            big_ind_arr = []
                            big_ind_arr.append(ind_arr)
                        last_wind = uplink_wind[-1][1]
                        if len(big_ind_arr) > 0:
                            D_1 = 0
                            D_65 = 0
                            D_97 = 0
                            D_113 = 0
                            D_121 = 0
                            D_125 = 0
                            temp = np.array(big_ind_arr).reshape(1, -1)
                            for o in temp[0]:
                                if o == 0:
                                    D_1 += 1
                                elif o == 64:
                                    D_65 += 1
                                elif o == 96:
                                    D_97 += 1
                                elif o == 112:
                                    D_113 += 1
                                elif o == 120:
                                    D_121 += 1
                                elif o == 124:
                                    D_125 += 1
                            uplink_wind[-1][2:len(uplink_wind[-1])] = (st.mode(temp[0]), D_1, D_65, D_97, D_113, D_121, D_125)
                    elif i - num_accum < last_wind:
                        uplink_wind[-1][1] = i + front_buf
                        last_wind = uplink_wind[-1][1]
                        if len(big_ind_arr) > 0:
                            D_1 = 0
                            D_65 = 0
                            D_97 = 0
                            D_113 = 0
                            D_121 = 0
                            D_125 = 0
                            temp = np.array(big_ind_arr).reshape(1, -1)
                            for o in temp[0]:
                                if o == 0:
                                    D_1 += 1
                                elif o == 64:
                                    D_65 += 1
                                elif o == 96:
                                    D_97 += 1
                                elif o == 112:
                                    D_113 += 1
                                elif o == 120:
                                    D_121 += 1
                                elif o == 124:
                                    D_125 += 1
                            uplink_wind[-1][2:len(uplink_wind[-1])] = (st.mode(temp[0]), D_1, D_65, D_97, D_113, D_121, D_125)
                t_n = 0
                t_p = 0
            if np.mod(c, num_accum) == 0:
                ind_arr = []
            else:
                ind_arr.append(ind)
            c += 1
            t_n = t_n + noise_floor
            t_p = t_p + fft_peak

            lst_ind = ind
        if len(uplink_wind) > 0:
            for i in range(0, len(uplink_wind)):
                uplink_wind[i][0:2] = np.multiply(uplink_wind[i][0:2], [win_jump, win_jump])
            temp = []
            for i in range(0, len(uplink_wind)):
                if ((uplink_wind[i][1] - uplink_wind[i][0]) / (windsz * upsampling_factor)) > 30:
                    temp.append(uplink_wind[i])
            uplink_wind = temp

        # if not self.nfl_flag:
            # self.chunk_num += 1

        # if not self.nfl_flag and self.chunk_num == 2:
            # # print(f"{max(peak_gain)}\t{self.nfl_flag}")
            # self.nfl = max(peak_gain) * 1.09
            # # self.nfl = max(peak_gain) * 6
            # # self.nfl = np.mean(peak_gain)
            # print(f"Your threshold is : {self.nfl}\n")
            # # self.nfl = max(peak_gain) * 1
            # self.nfl_flag = 1

        # self.lst_flag_right = self.cur_flag_right
        # self.lst_flag_left = self.cur_flag_left
        if len(uplink_wind) > 0:
            if uplink_wind[0][0] == win_jump:
                self.cur_flag_left = 1
            else:
                self.cur_flag_left = 0

            if (x_1_len - uplink_wind[-1][1]) < 0:
                self.cur_flag_right = 1
            else:
                self.cur_flag_right = 0
                self.chunk = x_1[uplink_wind[-1][0] : x_1_len]
        # print(f"Maximum of current peak gain is = {max(peak_gain)}")
        # print("big_ind_arr ", big_ind_arr)
        return uplink_wind
        
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
        
        if (ind < ((num_preamble)*window_size) ) :
            print("check_doppler_preamble edge case")
            shifts_small = np.arange(-ind, int(upsamp),1)
            shifts = np.arange(-ind, nshifts, int(upsamp))
        else:
            shifts_small = np.arange(-int(upsamp), int(upsamp),1)
            shifts = np.arange(-nshifts, nshifts, int(upsamp))
        out_preamb= 0
        possible_dop=0
        max_shift_arr = []
        max_bin_arr = []
        npt=32
        ind_old = ind
        max_val =  []
        for shift in shifts:
        # for i in [0]:
            ind_win = ((num_preamble)*window_size) + ind + shift
            win2 = np.concatenate(( np.conjugate(self.sym_2_css([0], N, SF, BW, FS)),self.sym_2_css([0], N, SF, BW, FS)))
            win1 =  Rx_Buffer[ind_win-window_size:ind_win+window_size]
            try: 
            # win2 = np.conjugate(self.sym_2_css([0], N, SF, BW, FS))
            # win1 =  Rx_Buffer[ind_win-window_size:ind_win]            
                doppler_fft = np.abs(np.fft.fft(win1 * win2) )
            except: 
                print("Failed INDICES:" , ind_win-window_size, ind_win+window_size)
                print("Ind: ", ind, "Queue len:", len(Rx_Buffer))
            max_val.append(max(doppler_fft))
        
        # plt.figure(5)
        # plt.plot(shifts, max_val)
        # plt.show()
        
        st_dx = shifts[np.argmax(max_val)]
        
        max_val = []
        argmax = []
        
        max_val2 = []
        
        sync_val = []
        
        doppler_acc_value = []
        doppler_acc_value2 = []
        
        sync_word_gain = []
        conj4_gain = []
        conj2_fft1_gain = []
        conj2_fft2_gain = [] 
        conj2add_gain = []
        conj2mult_gain = []
        doppler_dc_max = []
        doppler_dc_gain = []
        # plt.figure(15)
        # plt.plot(Rx_Buffer)
        # plt.axvline(ind + st_dx + min(shifts_small))
        # plt.axvline(ind + st_dx + max(shifts_small) + window_size, color='r')
        # plt.show()
        dop_argmax1 = []
        dop_argmax2 = []
        dop_argmax = []
        no_sync_doppler_fft_gain = []
        no_sync_val = []
        ''' Put a coarse-grained doppler shift here ''' 
        
        ind_win = ((num_preamble)*window_size) + ind + st_dx
        # test_dop = Rx_Buffer[ind_win-4*window_size:ind_win-2*window_size]\
        for i in range(0,10):
            test_dop = Rx_Buffer[ind_win-7*window_size+i:ind_win-1*window_size+i]
            window_time = int(len(test_dop)/2)
            test_dop = test_dop[::upsamp]
            est_doppler_slope = doppler_slope_offset2(test_dop, SF, BW, int(FS/upsamp), N, 0, 0, DC=False)
            est_doppler_slope = est_doppler_slope / (window_time/FS)
            # est_doppler_slope = 50
        print("Estimated doppler slope here: ", est_doppler_slope)
        # freq_list = fft_dx * (self.FS/(self.N*self.UPSAMP*npt))
        tt = np.arange(0, 3*window_size)/FS
        for shift in shifts_small:
        # for shift in np.arange(-10000, 10000):
        # for i in [0]:
            ind_win = ((num_preamble)*window_size) + ind + st_dx + shift
            # plt.figure(1)
            # plt.plot(self.sym_2_css([0], 2**(SF-1), (SF-1), int(BW*2), FS))
            # plt.show()
            
            # print(len(self.sym_2_css([0], 2**(SF+2), (SF+2), int(BW*2), FS)))
            # print(len(self.sym_2_css([0], 2**(SF), (SF), int(BW), FS)))
            
            # win2 = np.concatenate(( np.conjugate(self.sym_2_css([0], 2**(SF+2), (SF+2), int(BW*2), FS)),self.sym_2_css([0,0], 2**(SF+2), (SF+2), int(BW*2), FS)))
            win2 = np.concatenate(( np.conjugate(self.sym_2_css([0], 2**(SF), (SF), int(BW), FS)),self.sym_2_css([0,0], 2**(SF), (SF), int(BW), FS)))
            win1 =  Rx_Buffer[ind_win-window_size:ind_win+2*window_size] 
            # win1 =  np.concatenate(( Rx_Buffer[ind_win-window_size:ind_win+2*window_size] , np.zeros(len(Rx_Buffer[ind_win-window_size:ind_win+2*window_size])) )) 
            # win1 =  Rx_Buffer[ind_win:ind_win+window_size] 
            '''^This one DOES NOT cancel doppler'''
            # win2 = np.conjugate(self.sym_2_css([0], N, SF, BW, FS))
            # win1 =  Rx_Buffer[ind_win-window_size:ind_win]            
            no_sync_doppler_fft = np.abs(np.fft.fft(win1 * win2, len(win2)) )
            no_sync_val.append(max(no_sync_doppler_fft))
            no_sync_doppler_fft_gain.append(self.peak_gain_calculation(no_sync_doppler_fft))
            win1 = win1 * np.exp(1j * 2 * math.pi * -est_doppler_slope * (tt**2)/2)
            doppler_fft = np.abs(np.fft.fft(win1 * win2, len(win2)) )
            sync_val.append(max(doppler_fft))
            sync_fft = doppler_fft
            # win1 =  Rx_Buffer[ind_win-window_size:ind_win]
            # win2 =  Rx_Buffer[ind_win:ind_win+window_size]
            # win1 =  Rx_Buffer[ind_win-2*window_size:ind_win]
            # win2 =  Rx_Buffer[ind_win:ind_win+2*window_size]
            
            ''' This one cancels doppler '''
            win1 =  Rx_Buffer[ind_win-8*window_size:ind_win-4*window_size]
            win2 =  np.conjugate(Rx_Buffer[ind_win-4*window_size:ind_win])
            win3 =  Rx_Buffer[ind_win:ind_win+window_size]
            win4 = np.conjugate(Rx_Buffer[ind_win+window_size:ind_win+2*window_size])
            doppler_dc_fft = np.abs(np.fft.fft(win3*win4 ,int(len(win3))))
            doppler_dc_max.append(max(doppler_dc_fft))
            doppler_dc_gain.append(self.peak_gain_calculation(doppler_dc_fft))
            doppler_fft = np.abs(np.fft.fft(win1 * win2, int(len(win2))) ) # check if there is any gain in higher res FFT 
            # dop_max = max(doppler_fft)
            # dop_arg_max= np.argmax(doppler_fft)
            max_val.append(max(doppler_fft))
            # max_val.append(doppler_fft[0])
            argmax.append(np.argmax(doppler_fft))
            # win1n = Rx_Buffer[ind_win-8*window_size:ind_win-4*window_size]
            # win2n
            # win = Rx_Buffer[ind_win-4*window_size:ind_win]
            # win11 = win1n[:2*window_size]
            # win12 = np.conjugate(win1n[2*window_size:])
            # win21 = win2n[:2*window_size]
            # win22 = np.conjugate(win2n[2*window_size:])
            win11 = Rx_Buffer[ind_win-8*window_size:ind_win-6*window_size]
            win12 = np.conjugate(Rx_Buffer[ind_win-6*window_size:ind_win-4*window_size])
            win21 = Rx_Buffer[ind_win-4*window_size:ind_win-2*window_size]
            win22 = np.conjugate(Rx_Buffer[ind_win-2*window_size:ind_win])
                        
            # plt.figure(1)
            # plt.specgram(win11)
            # plt.figure(2)
            # plt.specgram(win12)
            # plt.figure(3)
            # plt.specgram(win21)
            # plt.figure(4)
            # plt.specgram(win22)
            
            # plt.show()
            
            # win11 = win1[:2*window_size]
            # win12 = np.conjugate(win1[2*window_size:])
            # win21 = win2[:2*window_size]
            # win22 = np.conjugate(win2[2*window_size:])
            
            fft_sz = int(len(win11)) * 1
            doppler_fft1 = np.fft.fft(win11 * win12, fft_sz)
            doppler_fft2 = np.fft.fft(win21 * win22, fft_sz)
            dop_argmax1.append(np.argmax(np.abs(doppler_fft1)))
            dop_argmax2.append(np.argmax(np.abs(doppler_fft2)))
            dop_argmax.append(np.argmax(doppler_fft))
            doppler_mult = np.abs(doppler_fft1 *doppler_fft2)
            doppler_add = np.abs(doppler_fft1 + doppler_fft2)
            # doppler_fft1 = np.abs(np.fft.fft(win11 * win12, fft_sz*16))
            # doppler_fft2 = np.abs(np.fft.fft(win21 * win22, fft_sz*16))
            # doppler_fft1 = np.abs(np.fft.fft(win11 * win21, fft_sz*16))
            # doppler_fft2 = np.abs(np.fft.fft(win12 * win22, fft_sz*16))
            # print("Peak gain at this point (mult two FFTs) : ",  self.peak_gain_calculation(doppler_mult))
            # print("Peak gain at this point (add two FFTs) : ",  self.peak_gain_calculation(doppler_add))
            doppler_acc_value.append(np.max(doppler_mult)) 
            doppler_acc_value2.append(np.max(doppler_add))
            # doppler_acc_value.append(doppler_mult[0]) 
            # doppler_acc_value2.append(doppler_add[0])
            dop_max1 = max(doppler_fft1)
            dop_max2 = max(doppler_fft2)
                 
            max_val2.append(dop_max1 * dop_max2)
            
            sync_word_gain.append(self.peak_gain_calculation (sync_fft) ) 
            conj4_gain.append(self.peak_gain_calculation (doppler_fft))
            conj2_fft1_gain.append(self.peak_gain_calculation (np.abs(doppler_fft1)))
            conj2_fft2_gain.append(self.peak_gain_calculation (np.abs(doppler_fft2)))
            conj2add_gain.append(self.peak_gain_calculation (doppler_add))
            conj2mult_gain.append(self.peak_gain_calculation (doppler_mult))
            
            # print("Max 4 chirps: ", dop_max, "Max 2 chirps: ", dop_max1, "Max 2 chirps: ", dop_max2)
            # print("Peak gains: ", peak_gain, peak_gain1, peak_gain2)
            '''
            Reference downchirp      
            
            '''            
            # win2 = np.conjugate(self.sym_2_css([0], N, SF, BW, FS))
            # win1 =  Rx_Buffer[ind_win-2*window_size:ind_win-window_size]            
            # doppler_fft2 = np.abs(np.fft.fft(win1 * win2, 32*len(win2)) )
            
            # max_val.append(max(doppler_fft) )
            # max_val.append(doppler_fft[0])
            # max_val.append(np.mean(doppler_fft[np.arange(-4,4)]))
        # abs_samp = np.argmax(max_val)
        abs_samp = np.argmax(sync_val)
        shift_dx = abs_samp        
        # plt.figure(1)
        # plt.plot(dop_argmax1)
        # plt.axvline(abs_samp, color='r')
        # plt.title('dop_argmax1')
        # # plt.show()
        
        # plt.figure(2)
        # plt.plot(dop_argmax2)
        # plt.axvline(abs_samp, color='r')
        # plt.title('dop_argmax2')
        
        
        # plt.figure(3)
        # plt.plot(dop_argmax)
        # plt.axvline(abs_samp, color='r')
        # plt.title('dop_argmax 4x4')
        # plt.show()
        
        if (False):
            # plt.figure(1)
            # plt.plot(max_val)
            # plt.xlabel('Samples Shift')
            # plt.axvline(abs_samp, color='r')
            # plt.ylabel('Max of FFT Dechirped at this Point')
            # plt.title('4 x conj[4]')
            
            # # plt.figure(2)
            # # plt.plot(max_val2)
            # # plt.axvline(abs_samp, color='r')
            # # plt.xlabel('Samples Shift')
            # # plt.ylabel('Max of FFT Dechirped at this Point')
            # # plt.title('( 2 x conj[2] ) * ( 2 x conj[2] )')
            
            # plt.figure(3)
            # plt.plot(sync_val)
            # plt.axvline(abs_samp, color='r')
            # plt.xlabel('Samples Shift')
            # plt.ylabel('Max of FFT Dechirped at this Point')
            # plt.title('2 DC + 1 UC With Doppler Slope Correction')
            # plt.figure(22)
            # plt.plot(no_sync_val)
            # plt.axvline(abs_samp, color='r')
            # plt.xlabel('Samples Shift')
            # plt.ylabel('Max of FFT Dechirped at this Point')
            # plt.title('2 DC + 1 UC without Doppler Slope Correction')
            
            # plt.figure(4)
            # plt.plot(doppler_acc_value)
            # plt.axvline(abs_samp, color='r')
            # plt.xlabel('Samples Shift')
            # plt.ylabel('Max of FFT Dechirped at this Point')
            # plt.title('Mult Two FFTs, THEN Take The Max')
            
            
            # plt.figure(5)
            # plt.plot(doppler_acc_value2)
            # plt.axvline(abs_samp, color='r')
            # plt.xlabel('Samples Shift')
            # plt.ylabel('Max of FFT Dechirped at this Point')
            # plt.title('Add Two FFTs, THEN Take The Max')
            
            plt.figure(66)
            plt.plot(no_sync_doppler_fft_gain)
            plt.axvline(abs_samp, color='r')
            plt.xlabel('Samples Shift')
            plt.ylabel('Gain (dB)')
            plt.title('Upchirp/Downchirp (without doppler slope suppression) GAIN ')
            
            
            plt.figure(6)
            plt.plot(sync_word_gain)
            plt.axvline(abs_samp, color='r')
            plt.xlabel('Samples Shift')
            plt.ylabel('Gain (dB)')
            plt.title('Upchirp/Downchirp GAIN ')
            
            plt.figure(7)
            plt.plot(conj4_gain)
            plt.axvline(abs_samp, color='r')
            plt.xlabel('Samples Shift')
            plt.ylabel('Gain (dB)')
            plt.title('4x4 Gain')
            
            plt.figure(8)
            plt.plot(conj2_fft1_gain)
            plt.axvline(abs_samp, color='r')
            plt.xlabel('Samples Shift')
            plt.ylabel('Gain (dB)')
            plt.title('2x2 First Half Gain')
            
            plt.figure(9)
            plt.plot(conj2_fft2_gain)
            plt.axvline(abs_samp, color='r')
            plt.xlabel('Samples Shift')
            plt.ylabel('Gain (dB)')
            plt.title('2x2 Second Half Gain')
            
            plt.figure(10)
            plt.plot(conj2add_gain)
            plt.axvline(abs_samp, color='r')
            plt.xlabel('Samples Shift')
            plt.ylabel('Gain (dB)')
            plt.title('2x2 + 2x2')
            
            plt.figure(11)
            plt.plot(conj2mult_gain)
            plt.axvline(abs_samp, color='r')
            plt.xlabel('Samples Shift')
            plt.ylabel('Gain (dB)')
            plt.title('2x2 * 2x2')
            
            # plt.figure(12)
            # plt.plot(np.add(conj2_fft1_gain , conj2_fft2_gain))
            # plt.axvline(abs_samp, color='r')
            # plt.xlabel('Samples Shift')
            # plt.ylabel('Gain (dB)')
            # plt.title('Directly add Peak Gains')
            # conj2mult_gain = []
            
            # plt.figure(12)
            # plt.plot(doppler_dc_gain)
            # plt.axvline(abs_samp, color='r')
            # plt.xlabel('Samples Shift')
            # plt.ylabel('Gain (dB)')
            # plt.title('doppler_dc_gain')
            # plt.figure(13)
            # plt.plot(doppler_dc_max)
            # plt.axvline(abs_samp, color='r')
            # plt.xlabel('Samples Shift')
            # plt.ylabel('Gain (dB)')
            # plt.title('doppler_dc_max')
           
            
            plt.show()
            # plt.figure(5)
        # plt.plot(shifts_small, max_val)
        # plt.show()
        
        # plt.figure(6)
        # plt.plot(shifts_small,argmax)
        # plt.show()
        
        
        # shift_small = shifts_small[np.argmax(max_val)]
        shift_small = shifts_small[shift_dx]
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
        

        # ind = 102397 #[REMOVE] hard-coded start valye
        
        # print("Current preamble shift: ", shifts[np.argmax(max_shift_arr)] )
        ind_win = ((num_preamble)*window_size) + ind #+ window_size
        # win2 = ( np.conjugate(self.sym_2_css([0], N, SF, BW, FS)))
        # win1 =  Rx_Buffer[ind_win-window_size:ind_win]
        win2 = ( np.conjugate(self.sym_2_css([0], N, SF, BW, FS)))
        # win1 =  Rx_Buffer[(ind_win-8*window_size):(ind_win-7*window_size)]
        win1 =  Rx_Buffer[ind:ind + window_size]
        npt = 128*4
        temp_wind_fft = abs(
                np.fft.fft(win1 * win2, n=npt*len(win2), axis=0))        
        max_bin = np.argmax(temp_wind_fft)
        
        # fft_dx = np.arange(-int(npt*len(win2)/2),int(npt*len(win2)/2))
        fft_dx = np.concatenate(( np.arange(0, int(npt*len(win2)/2)), np.arange(-int(npt*len(win2)/2),0)))
        max_bin = fft_dx[max_bin]
        # plt.figure(1)
        # plt.plot(win1 * win2)
        # plt.show()
        # if (max_bin < ((self.N*self.UPSAMP*npt)/2)):    
        # if (max_bin < ((self.N*self.UPSAMP*npt)/2)):    
            # max_bin = -((self.N*self.UPSAMP*npt/2) - max_bin)
        # freq_shift = self.FS * (max_bin/(self.N*self.UPSAMP*npt)
        freq_shift = max_bin * (self.FS/(self.N*self.UPSAMP*npt))
        freq_list = fft_dx * (self.FS/(self.N*self.UPSAMP*npt))
        # plt.figure(4)
        # plt.plot(fft_dx, freq_list)
        # plt.show()
        # freq_shift = 0
        # print(freq_shift, "freq shift preabmle", "Freq resolution: ", self.FS * 1/(self.N*self.UPSAMP*npt), "samp len:", len(win1))
        tt = np.arange(0,len(win1))/FS
        # win1=win1*np.exp(1j * 2 * math.pi * float(-freq_shift) * tt)    
        temp_wind_fft = abs(
                np.fft.fft(win1 * win2, n=npt*len(win2), axis=0))        
        # max_bin = np.argmax(temp_wind_fft)
        # print("INDEX OF THIS DOP:", ind)
        # plt.figure(2)
        # plt.plot(temp_wind_fft[np.arange(-40,40)])
         
        # ffine= FFO_corr_sym(win1, SF, BW, FS, N, 0)
        ffine = 0 #[TODO] replace this with the appropraite values
        # plt.show()    
        # possible_dop=int(freq_shift)
        possible_dop=float(freq_shift+ffine)
        # possible_dop = 0 #[TODO] remove
        # print("Possible frequency shift:", freq_shift)
        out_preamb=ind 
         
        # print("Preamble doppler time: ", time.time() - t_start)
        return out_preamb, possible_dop         

    def peak_gain_calculation(self,abs_fft):
        '''
            abs_fft is absolute value of dechirped FFT 
        '''
        
        fft_peak = np.max(abs_fft)
        # fft_peak = abs_fft[0]
        max_arg = np.argmax(abs_fft) 
        fft_dx = np.concatenate(( np.arange(0, max_arg),np.arange(max_arg+1, len(abs_fft)))) 
        # fft_dx = np.arange(1, len(abs_fft))
        noise_obs = abs_fft[fft_dx]
        noise_floor = np.mean(noise_obs)            
        
        # fft_peak = np.max(doppler_fft)            
        peak_gain = 10 * math.log10(fft_peak/noise_floor) 
        # peak_gain = abs_fft[0]
        return peak_gain
        
    def demod_preamble(self, queue, ref_preamble, window_size, FS,freq_shift1,FD_COARSE,FD_FINE,N,SF,BW):
        '''
        Get dopppler from the preamble, and decode packet length. 
        '''   
        # print(window_size, int(FS/BW) * N)
        # print(len(queue), window_size*10)
        # give a few samples prior to the estimated preamble location to account for fine sampling offset     
        ind = 0 
        tt = np.arange(0, len(queue))/FS 
        num_preamble = 8
        idx_shift= FFO_corr(queue*np.exp(1j * 2 * math.pi * float(freq_shift1) * tt),ind, SF, BW, FS, N, num_preamble)
        idx_shift = 0 # TODO [REMOVE] 
        queue = queue[idx_shift:]
        # print(freq_shift1, "Shift at preamble.")
        # print("FFO shift:", idx_shift)
        
        demod_status = False     
        
        tp = np.linspace(0, (len(ref_preamble)+2*window_size)/FS, len(ref_preamble)+2*window_size)
        
        
        freq_shift = freq_shift1    
        f_old = freq_shift
        ts = np.linspace(0, window_size/FS, window_size)
        # tt = np.linspace(0, window_size/FS, window_size)
        tt = np.arange(0, 2*window_size)/FS
        # check to make sure we can demodulate the preamble. If not, something went wrong with our sync and we need to try another pt 
        
        
        self.curr_sym = 0 
        self.prev_sym = 0
        # corr_dopp = np.concatenate((self.UPCHIRP, queue[:window_size]* np.exp(1j * 2 * math.pi * freq_shift * ts)))
        # freq_shift2 = doppler_slope_offset(corr_dopp, self.SF, self.BW, self.FS, self.N, self.prev_sym, self.curr_sym)     #[TESTING]               # f_shift_new=0 #[TODO] Remove.
        # self.prev_sym = self.curr_sym 
        # freq_shift = freq_shift + freq_shift2
        
        
        # self.curr_sym = 0 
        # self.prev_sym = 0
        # corr_dopp = np.concatenate((self.UPCHIRP, queue[window_size*(0):window_size*(1)]* np.exp(1j * 2 * math.pi * freq_shift * ts)))
        # # freq_shift2 = doppler_slope_offset(corr_dopp, self.SF, self.BW, self.FS, self.N, self.prev_sym, self.curr_sym)     #[TESTING]               # f_shift_new=0 #[TODO] Remove.
        # # self.prev_sym = self.curr_sym 
        # # freq_shift = freq_shift + freq_shift2P, queue[window_size*(1):window_size*(1+1)]* np.exp(1j * 2 * math.pi * freq_shift * ts)))
        # freq_shift2 = doppler_slope_offset(corr_dopp, self.SF, self.BW, self.FS, self.N, self.prev_sym, self.curr_sym)     #[TESTING]               # f_shift_new=0 #[TODO] Remove.
        # # self.prev_sym = self.curr_sym 
        # print("Small initial frequency shift: " , freq_shift2)
        # freq_shift = freq_shift + freq_shift2
        
        # print("Initial doppler offset ocmpensation/")
        
        
        # print("Demodulated first symbol: ", self.symbol_demod_sig(queue[window_size*(1):window_size*(1+1)]* np.exp(1j * 2 * math.pi * freq_shift * ts)) )
        doppler_sample = np.ones(window_size)
        self.prev_sym = 0
        d_slope_avg = []
        for symdx in range(1,8):
            # print(self.symbol_demod_sig(self.PREV_QUEUE[self.WINDOW_SIZE*(symdx):self.WINDOW_SIZE*(symdx+1)]))
            # doppler = (np.arange(0, window_size) / FS ) * 5
            # sym_win = queue[window_size*(symdx):window_size*(symdx+1)] * np.exp(1j * 2 * math.pi * doppler * ts)
            # print(freq_shift, "Shift at sym.")
            sym_win = queue[window_size*(symdx):window_size*(symdx+1)] * np.exp(1j * 2 * math.pi * freq_shift * ts)
            # sym_win = queue[window_size*(symdx):window_size*(symdx+1)]
            # sym_win = queue[window_size*(symdx):window_size*(symdx+1)] * np.exp(1j * 2 * math.pi * freq_shift * ts) * np.conjugate(doppler_sample)
            dem_sym = self.symbol_demod_sig(sym_win) 
            # freq_shift2 = FFO_corr_sym(sym_win, self.SF, self.BW, self.FS, self.N, 0)
            
            self.curr_sym = dem_sym 
            freq_shift2 = -doppler_slope_offset2(queue[window_size*(symdx-1):window_size*(symdx+1)]* np.exp(1j * 2 * math.pi * freq_shift * tt), self.SF, self.BW, self.FS, self.N, self.prev_sym, self.curr_sym)     #[TESTING]               # f_shift_new=0 #[TODO] Remove.
            self.prev_sym = self.curr_sym 
            d_slope_avg.append(freq_shift2) 
            freq_shift = freq_shift + freq_shift2
            # print("New f shift: " , freq_shift2, "Total f shift: ",freq_shift)
            # print("est. doppler: ",freq_shift2 / tt[window_size])
            if (dem_sym != 0):
                print("Sync or doppler correction is incorrect.",self.symbol_demod_sig(sym_win))
                print("Sym dx: ", symdx)
                demod_status = False
                return demod_status, 0,0,0  
        # print(d_slope_avg)
        self.doppler_slope = np.mean(d_slope_avg)
        # print("AVERAGE DOPPLER SLOPE : ", self.doppler_slope / tt[window_size],tt[window_size])
        # print("Doppler slope per sym: ", self.doppler_slope) 
        
        freq_shift = freq_shift + self.doppler_slope
        for downchirp in range(0,2): 
            # detect the downchirps (sync word) 
            # sym_win = queue[window_size*(downchirp+8):window_size*(downchirp+9)] * np.exp(1j * 2 * math.pi * freq_shift * ts)
            # sym_win = queue[window_size*8 + (window_size*downchirp):window_size*(downchirp+8)+(window_size*2*downchirp)] * np.exp(1j * 2 * math.pi * freq_shift * ts)
            st_dx = window_size*8 + downchirp*window_size 
            end_dx = st_dx+window_size
            print(st_dx, end_dx)
            sym_win = queue[st_dx:end_dx] * np.exp(1j * 2 * math.pi * freq_shift * ts)
            
            # plt.figure(1)
            # plt.specgram(sym_win)
            # plt.show()
            
            # dem_sym = self.symbol_demod_sig_dc(sym_win)
            dem_sym = self.symbol_demod_sig(np.conjugate(sym_win))            
            # print("Preamble symbol: ", dem_sym)
            # freq_shift2 = FFO_corr_sym_conj(sym_win, self.SF, self.BW, self.FS, self.N, 0)
            # freq_shift2 = 0
            # freq_shift2 = -FFO_corr_sym(np.conjugate(sym_win), self.SF, self.BW, self.FS, self.N, 0)
            # freq_shift2 = FFO_corr_sym(np.conjugate(sym_win), self.SF, self.BW, self.FS, self.N, 0) # [TODO] Remove   
            # freq_shift = freq_shift + freq_shift2
            freq_shift = freq_shift + self.doppler_slope         
            # freq_shift = 0 
            # print("New f shift [coj]: " , freq_shift2, dem_sym, "Total f shift: ",freq_shift)
            
            if (dem_sym != 0):
                print("Sync or doppler correction is incorrect. [DC]",self.symbol_demod_sig(sym_win))
                print("Sym dx: ", downchirp, dem_sym)
                demod_status = False
                return demod_status, 0,0,0  
            # else:
                # print("Passed.")       


                
        # t  = np.linspace(0, len(queue)/FS, len(queue))
        # queue = queue * np.exp(1j * 2 * math.pi * freq_shift * t)   
        self.OUTPUT_PKT = []
        end_preamb_dx = 2*window_size+len(ref_preamble)
        
        # freq_shift = freq_shift + freq_shift2        
        pkt_len_2 = self.symbol_demod_sig(queue[end_preamb_dx:end_preamb_dx+window_size]* np.exp(1j * 2 * math.pi * freq_shift * tt[:window_size]) )
        # print("first packet length: ", pkt_len_2)
               
        # self.curr_sym = pkt_len_2
        # self.prev_sym = 0
        # # dop_win = np.conjugate(queue[end_preamb_dx-window_size:end_preamb_dx]* np.exp(1j * 2 * math.pi * freq_shift * ts)  )
        # dop_win = queue[end_preamb_dx-window_size:end_preamb_dx]* np.exp(1j * 2 * math.pi * freq_shift * ts) 
        # dop_win =  np.concatenate(( dop_win, queue[end_preamb_dx:end_preamb_dx+window_size]* np.exp(1j * 2 * math.pi * freq_shift * ts )))

        
        # freq_shift2 = doppler_slope_offset(dop_win, self.SF, self.BW, self.FS, self.N, self.prev_sym, self.curr_sym,DC=True) 
        # freq_shift = freq_shift + 
        freq_shift = freq_shift + self.doppler_slope 
        self.prev_sym = self.curr_sym
        
        # freq_shift2 = FFO_corr_sym(queue[end_preamb_dx:end_preamb_dx+window_size]* np.exp(1j * 2 * math.pi * freq_shift * t[:window_size]), self.SF, self.BW, self.FS, self.N, pkt_len_2)
        # print("-New f shift p1: " , freq_shift2)
        # freq_shift = freq_shift+ freq_shift2
            
        pkt_len_1 = self.symbol_demod_sig(queue[end_preamb_dx+window_size:end_preamb_dx+2*window_size]* np.exp(1j * 2 * math.pi * freq_shift * tt[:window_size]))
        
        self.curr_sym = pkt_len_1
        
        # dop_win =  np.concatenate(( dop_win, queue[end_preamb_dx:end_preamb_dx+window_size]* np.exp(1j * 2 * math.pi * freq_shift * ts] ))          
        # freq_shift2 = doppler_slope_offset(queue[end_preamb_dx:end_preamb_dx+2*window_size]* np.exp(1j * 2 * math.pi * freq_shift * tt), self.SF, self.BW, self.FS, self.N, self.prev_sym, self.curr_sym) 
        # freq_shift = freq_shift + freq_shift2
        freq_shift = freq_shift + self.doppler_slope
        
              
        self.prev_sym = pkt_len_1
        # freq_shift2 = FFO_corr_sym(queue[end_preamb_dx+window_size:end_preamb_dx+2*window_size]* np.exp(1j * 2 * math.pi * freq_shift * t[:window_size]), self.SF, self.BW, self.FS, self.N, pkt_len_1)
        # print("-New f shift p2: " , freq_shift2)
        # freq_shift = freq_shift+ freq_shift2       
                        
        #[TODO] : What if other candidate IDXs lie within this, and this sync is just a bit off? Need to not remove the pts then
        packet_len =  pkt_len_2 + N*pkt_len_1
        if (packet_len > 1000):
            print(" ===== found pkt length too long ", packet_len, pkt_len_2, pkt_len_1)
            pass
        else: 
            demod_status = True    
        # print("Packet length:", packet_len)    
        # print(" === END OF PREAMBLE DEMOD === ")
        # print("Freq shift at preamble : ", freq_shift - f_old, f_old)
        
        # print("Selected IDX: ",idx_shift)
        
        return demod_status, packet_len, freq_shift, idx_shift        

    
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
            # print("Ground Truth:",self.GND_TRUTH_PKT)
            # print("Output:      ", self.OUTPUT_PKT)
            if (len(self.OUTPUT_PKT) == len(self.GND_TRUTH_PKT)):
                nwrong = np.count_nonzero(np.abs(np.subtract(self.OUTPUT_PKT , self.GND_TRUTH_PKT)))
                print("num wrong: ",nwrong)    
                if (nwrong > 0):
                    print('\007')                
            else: 
                print("Packet length does not agree." )
                print('\007')  
            print("TOTAL PACKET COUNT: ", self.TOTAL_PKT_CNT)
            print("SNR For this packet:",self.PKT_SNR)
        else: 
            print("Output:      ", self.OUTPUT_PKT)
            with open(self.OUTFILE+'_'+str(self.TOTAL_PKT_CNT) + '.pkl', 'wb') as f: 
                pickle.dump([self.GND_TRUTH_PKT,  self.OUTPUT_PKT, self.TOTAL_PKT_CNT,self.PKT_SNR, self.CURR_POSIX_TIME],f)
            
