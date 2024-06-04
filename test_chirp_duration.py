'''
Code for testing CSS chirp generation
Last modified 6/4/2024 

'''

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

def sym_2_css(symbol, N, SF, BW, Fs, upsamp): 
    '''
    sym_2_css :: Modules a symbol between (0,N-1) to a chirp  
    
    N = 2^SF 
    SF 
    BW = Bandwidth sweep of the chirp -BW/2 to BW/2 
    Fs = Sampling rate
    '''
    sym_out = []
    
    # Fs/Bw is upsamp...
    # spsym = int(Fs/BW*N) # symbols defined by their freq offset at the start 
    # spsym = Ts
    spsym = N*upsamp
    print("usamp vs samp:",N*upsamp, int(Fs/BW*N))
    # spsym = N*10
    
    T = spsym/Fs
    k = BW/(T) 
    # k = BW/spsym 
    f_start = -BW/2 
    pi = math.pi 

    for sym in symbol: 
        st_offset = int(sym)
        # t_arr = np.arange(0, spsym*(N-st_offset)/N)/Fs
        # t_arr = np.arange(0, spsym*(N-st_offset)/N)/Fs
        t_arr = np.arange(0, int(spsym*(N-st_offset)/N))/Fs
        # cs = np.exp(1j * 2 * pi *(t_arr*(f_start + k*T*st_offset/N + 0.5*k*t_arr)))   # don't forget residuals
        print(BW*(st_offset/N),"off")
        cs = np.exp(1j * 2 * pi *(t_arr*(f_start + BW*(st_offset/N) + 0.5*k*t_arr)))
        if (len(t_arr) >0):
            # ang_start = np.angle(cs[len(t_arr)-1]) # we need to keep track of the phase 
            ang_start = np.angle(cs[-1])
        else: 
            ang_start = 0
                
        #t_arr = np.arange(0, spsym*st_offset/N - 1)/Fs # [TODO] fix this ... 
        # t_arr = np.arange(0, spsym*st_offset/N )/Fs 
        # t_arr = np.arange(0, spsym*st_offset/N )/Fs 
        t_arr = np.arange(0, int(spsym*st_offset/N ))/Fs 
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
def symbol_demod_sig(sig_tmp, UPCHIRP, UPSAMP,N):
    '''
    Demodulates a CSS symbol and returns the frequency bin 
    at which the symbol appears.
    '''        
    
    trans_upchirp = np.conjugate(UPCHIRP)
    
    # sig_tmp = scipy.signal.decimate(sig_tmp, UPSAMP)
    # trans_upchirp = scipy.signal.decimate(trans_upchirp, UPSAMP)
    sig_tmp = sig_tmp[::UPSAMP]
    # trans_upchirp = scipy.signal.decimate(trans_upchirp, UPSAMP)
    trans_upchirp = trans_upchirp[::UPSAMP]
    
    # trans_upchirp = trans_upchirp[::10]
    dechirped = sig_tmp * trans_upchirp
            
    dechirped = np.squeeze(dechirped)
    print(len(dechirped))
               
    
    # dechirped = dechirped[::10]

    # data_fft = abs(np.fft.fft(dechirped)).transpose()
    # data_fft = abs(np.fft.fft(dechirped))
    
    # print(len(data_fft), "dat fft")
    dechirped = abs(np.fft.fft(dechirped))
    # dechirped = np.concatenate((data_fft[:int(N/2)], \
         # data_fft[-int(N/2):]))
    # dechirped = data_fft
    print(len(dechirped))
    freq_bin = np.argmax(dechirped)    
        
    return freq_bin        

def create_upchirp():
    '''
    Create an upchirp which is used to demodulate the symbols 
    '''
    return sym_2_css([0], self.N, self.SF, self.BW, self.FS)
        
if __name__ == "__main__": 
    Fs=20000
    # Ts = 1280 # whatever 20KHz symbol period was... nSamples for a symbol Ts/Fs = symbol period 
    SF = 9 
    N = 2**SF
    BW = 10000
    upsamp = 1
    uchirp = sym_2_css([0], N, SF, BW, Fs,upsamp)
    
    # cc=sym_2_css([0], N, SF, BW, Fs,20)
    # cc = cc[::20]
    # cc2 = sym_2_css([0], N, SF, BW, Fs,upsamp)
    # cc2 = cc2[::10]
    # plt.figure(1)
    # plt.plot(cc)
    # plt.plot(cc2)
    # plt.show()
    # sym_tmp = sym_2_css([int(N/2)], N, SF, BW, Fs,upsamp) 
    sym_tmp = sym_2_css([64], N, SF, BW, Fs,upsamp) 
    
    demod_sym = symbol_demod_sig(sym_tmp, uchirp, upsamp,N)
    print("ANS: ", demod_sym)
    plt.figure(1)
    plt.specgram(uchirp)
    plt.show()
    # # sym_tmp = scipy.signal.decimate(sym_tmp, upsamp)
    # sym_tmp=sym_tmp[::upsamp]
    # # uchirp = scipy.signal.decimate(uchirp, upsamp)
    # uchirp = uchirp[::upsamp]
    # plt.plot(np.fft.fft(np.conjugate(uchirp) * sym_tmp) )
    # plt.show()