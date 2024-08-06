import math
import multiprocessing
from multiprocessing import Manager
import socket
import cmath
import time
from threading import Thread
from queue import Queue
import numpy as np
from utils import *
#from css_demod import CssDemod
from css_demod_2 import CssDemod
from FFO_corr import FFO_corr_UPSAMP 
import pickle
import zmq
from stft_v2 import stft_v2
from matplotlib import pyplot as plt

def sym_2_css(symbol, N, SF, BW, Fs): 
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

        t_arr = np.arange(0, spsym*st_offset/N )/Fs 
        ce = np.exp(1j*(ang_start + 2*pi*(t_arr*(f_start +0.5*k*t_arr)) ))
        out_sig = np.concatenate((cs, ce))
        sym_out=np.append(sym_out,out_sig)
       
    return sym_out

def FFO_corr_sym(Data_stack, SF, BW, Fs, N, sym):
    '''
    Determines the fine frequency offset from the received samples in Datastack to the 
    ideal generated symbol by dechirping with an upsampled FFT with resolution Fs/M and
    identifying a frequency offset. 
    '''
    upsampling_factor = int(Fs/BW)
    # print(Fs)
    # print(upsampling_factor)
    # DC = np.conjugate(sym_2_css([0], N, SF, BW, Fs))
    # true_sym = sym_2_css([sym], N, SF, BW, Fs)
    true_sym = np.conjugate(sym_2_css([sym], N, SF, BW, Fs))
    # npt = 16
    npt = 128
    # npt = 4
    # npt = 64
    
    # [TODO] Solve for base 2 number that provides the appropriate frequency resolution
    
    inn = []
    bin_vals = []
    freq_offsets = [] 
    # for i in np.arange(-5,5):
    for i in [0]:
        # freq_mean = [] 
        # for pp in range(num_preamble-1):
        data_wind = np.roll(Data_stack, i) # simulate different starting frequencies 
        data_wind = data_wind * np.exp(1j * -np.angle(data_wind[0]))
        
        # plt.figure(7)
        # plt.plot(np.angle(data_wind))
        # plt.show()
        # dech = np.abs(np.fft.fft(np.concatenate([(data_wind * true_sym))))
        # dech1 = np.abs(np.fft.fft(np.concatenate([(data_wind * true_sym)]), npt*len(data_wind)))
        dech1 = np.abs(np.fft.fft(np.concatenate([(data_wind * true_sym)]), npt*N*upsampling_factor))
        # dech1 = np.concatenate([dech1[:npt], dech1[-npt:]])
        dech1 = np.concatenate([dech1[-npt:],dech1[:npt]])
        # plt.figure(2)
        # plt.plot(dech1)
        # plt.show()
        
        bin_offset = np.argmax(dech1)
        bin_val = np.max(dech1) 
        bin_vals.append(bin_val) 
        if (bin_offset > npt):
            # print("Positive offset.")
            bin_offset =  bin_offset - npt 
        else:
            bin_offset = -(npt - bin_offset)
        freq_offset = (Fs) * bin_offset/(N*upsampling_factor*npt) #* 10
        # print("Resolution: ", (Fs) * 1/(N*upsampling_factor*npt))
        freq_offsets.append(freq_offset) 
        # print(bin_offset, "OFFSET")
    freq_offset = freq_offsets[np.argmax(bin_vals)]
    print("Tentative frequency offset: ", freq_offset)
        
    return -freq_offset     
def doppler_slope_offset(Data_stack, SF, BW, FS, N, prev_sym, curr_sym, DC=False):
    '''
    Determine doppler offset between two symbols for continuous doppler correction.
    '''
    us_factor = int(FS/BW)
    sym_size = us_factor * N
    fft_us = 16 # upsample the fft by this factor 
    nfft = fft_us * sym_size 
    fft_idx = np.arange(-int(nfft/2),int(nfft/2))
    # fft_idx = np.arange(0, int(nfft))
    
    # if (prev_sym > curr_sym): 
        # fft_idx = -np.arange(0, int(nfft/2)-1)
    # else:
    # fft_idx = np.arange(0, int(nfft/2)+1)

    if (DC):    # if one of the chirps is a downchirp (such as sync word)         
        plt.figure(6)
        plt.specgram(Data_stack[:sym_size])
        
        s1 = Data_stack[:sym_size] *  sym_2_css([0],N,SF,BW,FS) # must be flipped so that the doppler does not accumulate and cause spreading
        s2 = np.conjugate( Data_stack[sym_size:] * np.conjugate(sym_2_css([0],N,SF,BW,FS)) )
                        
        plt.figure(7)
        plt.plot(s1)
        # plt.specgram(s1)
        plt.figure(8)
        plt.plot(s2)
        
        plt.figure(3)
        plt.plot(np.fft.fft(s1 * s2))
                
        doppler_slope = s1 *s2 
        s1 = Data_stack[:sym_size] * (-doppler_slope) 
        s2 = Data_stack[sym_size:] * (-doppler_slope) 
        
        # plt.figure(8)
        # plt.specgram(freq_correction)
        plt.show()
        
        exp_idx = ( ( curr_sym-prev_sym ) * (fft_us))       
        exp_idx = 0
    else: 
        s1 = np.conjugate(Data_stack[:sym_size])
        s2 = Data_stack[sym_size:]
        # get   the index of perfect demodulation without doppler shift
        # expected 'tone' of dechirping one symbol with another (non-DC usually)
        exp_idx = ( curr_sym-prev_sym ) % (N*us_factor) * (fft_us) #* N
        exp_idx = exp_idx # FFT offset may be two-sided
    # print(( (curr_sym * us_factor) - (prev_sym * us_factor) ))
    # print(exp_idx)
    # exp_idx = fft_idx[exp_idx]
    # print("Uncorrected IDX" ,exp_idx, "Prev sym: ", prev_sym, "Curr sym: ", curr_sym)
    dechirped = np.abs(np.fft.fft(s1 * s2, nfft))
    
    # print(np.argmax(dechirped[np.arange(0, int(nfft/2))]), np.argmax(np.flip(dechirped[-np.arange(int(nfft/2), int(nfft)+1) ])))
    
    # plt.figure(1)
    # plt.plot(np.concatenate((np.arange(0, int(nfft/2)),  np.arange(int(nfft/2), int(nfft)))) ) 
    # plt.show()
    # pos_max = np.argmax(dechirped[np.arange(0, int(nfft/2))])
    # neg_max = np.argmax(np.flip(dechirped[-np.arange(int(nfft/2), int(nfft)+1) ]))
    
    # if (pos_max != exp_idx):
        # plt.figure(1)
        # plt.plot(dechirped)
        # plt.plot(fft_idx,dechirped)
        # plt.show()        
        
    # dechirped = dechirped[fft_idx]
        
    freq_dx = fft_idx * (FS / (nfft))
    # plt.figure(1)
    # plt.plot(freq_dx, dechirped)
    # plt.plot(fft_idx,dechirped)
    # plt.show()
    # doppler_offset = fft_idx[np.argmax(dechirped)]
    doppler_offset = np.argmax(dechirped)
    
    sym_dx = doppler_offset - (doppler_offset % (fft_us*us_factor))
    guess_sym = (sym_dx - prev_sym) % (N) 
    print("Guessed demodulated sym: ", guess_sym, "True sym: ", curr_sym)
    # doppler_offset = pos_max
    # doppler_offset = np.argmax(dechirped)
    # print("Doppler IDX Offset: ", doppler_offset, "Len: ", len(dechirped))
    if (prev_sym > curr_sym):
        frequency_offset = (exp_idx-doppler_offset) * (FS / (nfft))
    else: 
        frequency_offset = -(exp_idx-doppler_offset) * (FS / (nfft))
    print("Doppler IDX Offset: ", doppler_offset,"Prev sym: ", prev_sym, "Curr sym: ", curr_sym,"Current frequency offset: ", frequency_offset) 
    
    return frequency_offset
FS =200000
# FS =5000
# SLOPE_TIME = 5*60*FS # number of samples, minutes * seconds * samples/second 
# doppler_slope = 100 # Rate of change of doppler, Hz/second
doppler_slope = 50 # Rate of change of doppler, Hz/second
doppler_start = 0
# starting_frequency = 0
# starting_frequency = 0


SF =7 
N = 2**SF 
BW = 2500 
# BW = 20000
upsampling_factor = int(FS/BW)
# FS
symbol_size = len(sym_2_css([56],N,SF,BW,FS))*2

sampling_offset = 0
# sampling_offset = 0
time_vec = (np.arange(sampling_offset, symbol_size+sampling_offset) / FS )
doppler =  (((np.arange(0, symbol_size) / FS ) * doppler_slope) + doppler_start) 
start_freq = doppler[0]
end_freq = doppler[-1]
print("Overall ending frequency: ", end_freq)
# doppler_out = np.exp(1j*( 2*math.pi*(doppler *time_vec) ))
doppler_out = np.exp(1j * 2*math.pi *(doppler_start*time_vec+ doppler_slope * ((time_vec)**2)/2)) # like classic chirp

plt.figure(9)
plt.plot(doppler_out)
# plt.plot(np.angle(doppler_out2))
# plt.show()
# dop_ang=np.angle(doppler_out)

# sampling_offset = symbol_size*1
# time_vec = (np.arange(sampling_offset, symbol_size+sampling_offset) / FS )
# time_vec0 = (np.arange(0, symbol_size+0) / FS )
# doppler =  (((np.arange(0, symbol_size) / FS ) * doppler_slope) + doppler_start) 
# start_freq = doppler[0]
# end_freq = doppler[-1]
# print("Overall ending frequency: ", end_freq)
# doppler_outP = np.exp(1j*( 2*math.pi*(doppler *time_vec) ))
# doppler_outx = np.exp(1j*( 2*math.pi*(doppler *time_vec) )) * np.exp(1j*-dop_ang*time_vec0**2)
# # plt.figure(9)
# # plt.plot(np.angle(doppler_outP))
# plt.plot(np.angle(doppler_outx))
# plt.plot(np.angle(doppler_out2))
# plt.show()

''' Frequency change sanity check ''' 
Freq = (np.angle(doppler_out[1:] * np.conjugate(doppler_out[0:-1])));
Freq = Freq / (2*math.pi) * (FS);
# plt.figure(1)
# plt.plot(time_vec[0:-1],Freq)
# plt.show()
# plt.figure(9)
# plt.plot(np.round(doppler_out2 - doppler_out))
# plt.figure(8)
# plt.specgram(doppler_out)
# plt.show()
# print("Mean angle diff: ",np.mean(np.diff(np.angle(doppler_out))))
# symbol_train = sym_2_css([1,1],N,SF,BW,FS)
first_sym = 100
end_sym = 20


for first_sym in range(0, N-1):
    for end_sym in range(0,N-1):
# for first_sym in range(0,1):
    # for end_sym in range(0,N-1):
        symbol_train = sym_2_css([first_sym,end_sym],N,SF,BW,FS)
        
        # noise = np.random.normal(0,.5,len(symbol_train))
        
        # plt.figure(1)
        # plt.plot(symbol_train)
        # SNR = 10*np.log10(np.mean(np.abs(symbol_train)) / np.mean(np.abs(noise)))
        # print("Current SNR: ", SNR)
        # symbol_train = symbol_train + noise 
        
        
        # plt.figure(2)
        # plt.plot(symbol_train)
        
        # plt.show()
        
        # symbol_train = np.concatenate(( np.conjugate(sym_2_css([0],N,SF,BW,FS)),sym_2_css([0],N,SF,BW,FS) )) 
        # symbol_train = np.concatenate(( np.conjugate(sym_2_css([0],N,SF,BW,FS)),sym_2_css([0],N,SF,BW,FS) )) 
        # d_off = doppler_slope_offset(symbol_train* doppler_out, SF, BW, FS, N, first_sym, end_sym, DC=False)
        d_off = doppler_slope_offset(symbol_train[::upsampling_factor]* doppler_out[::upsampling_factor], SF, BW, FS/upsampling_factor, N, first_sym, end_sym, DC=False)

# d_off = doppler_slope_offset(symbol_train* doppler_out, SF, BW, FS, N, 30, 30)
# print("doppler Offset: ", d_off)



# DC = np.conjugate(symbol_train)
# npt = 16
# # fft_dx = np.arange(-int(len(DC)*npt/2), int(len(DC)*npt/2))
# # fft_dx = np.concatenate( (  np.arange(-int(len(DC)*npt/2), int(len(DC)*npt) ), np.arange(0,int(len(DC)*npt) )))


# # fft_dx = np.concatenate( (  np.arange(-int(len(DC)*npt/2), int(len(DC)*npt) ), np.arange(0,int(len(DC)*npt) )))
# # fft_dx = np.arange(0, len(DC)*npt)

# # doppler_fft = np.abs(np.fft.fft(doppler_out*symbol_train, npt*len(DC)))
# # doppler_fft = np.abs(np.fft.fft(symbol_train*doppler_out*DC, npt*len(DC)))
# # doppler_fft = np.abs(np.fft.fft(doppler_out, npt*len(DC)))
# # unflipped = symbol_train[int(symbol_size/2):]
# unflipped = symbol_train[0:int(symbol_size/2)] * doppler_out[0:int(symbol_size/2)] 
# # unflipped = np.arange(0,10)
# # flipped_pt2 = symbol_train[int(symbol_size/2):]

# doppler_sym = np.conjugate(symbol_train[int(symbol_size/2):] * doppler_out[int(symbol_size/2):])
# # flipped_pt2 = np.conjugate(doppler_sym[::-1])
# flipped_pt2 = doppler_sym#[::-1]
# # flipped_pt2 = np.zeros((len(unflipped)))
# # print(unflipped.shape)
# # print(flipped_pt2.shape)
# # freq_arr = []
# # for pt in range(0, len(unflipped)):
    # # flipped_pt2[pt] = unflipped[len(unflipped) - pt -1]
    # # freq_arr.append(len(unflipped) - pt -1)

# # plt.figure(8)
# # plt.plot(np.arange(0, len(unflipped)), freq_arr)
# # print(freq_arr)
# dout = flipped_pt2 * unflipped
# Freq = (np.angle(dout[1:] * np.conjugate(dout[0:-1])));
# # Freq = (np.angle(flipped_pt2[1:] * np.conjugate(flipped_pt2[0:-1])));
# Freq = Freq / (2*math.pi) * (FS);
# plt.figure(1)
# plt.plot(time_vec[0:int(symbol_size/2)-1],Freq)
# plt.show()

# # flipped_pt2 = np.flip(unflipped)
# # print(flipped_pt2.shape)
# # rev_ind = np.arange(0,-int(symbol_size/2),1)
# # print(rev_ind)
# # flipped_pt2 = symbol_train[rev_ind]
# plt.figure(1)
# plt.specgram(unflipped)
# plt.figure(2)
# plt.specgram(flipped_pt2)
# plt.show()
# # unflipped = symbol_train[:int(symbol_size/2)]
# doppler_fft = np.abs(np.fft.fft(flipped_pt2 * unflipped, npt*len(DC)))  
# # print(len(DC)*npt)
# fft_dx = np.arange(-int(len(doppler_fft)/2), int(len(doppler_fft)/2))
# # fft_dx = np.arange(-int(len(DC)*npt/2), int(len(DC)*npt/2))
# # print(len(fft_dx), len(doppler_fft))
# # plt.figure(7)
# # plt.plot(fft_dx)
# # plt.show()
# freq_values = fft_dx *((FS)/(npt*symbol_size)) 

# # doppler_fft = doppler_fft[fft_dx]
# doppler_fft = doppler_fft[fft_dx]
# freq_max = np.argmax(doppler_fft)
# print(freq_max)
# freq_max = freq_values[freq_max]
# print("max frqu: ", freq_max)
# plt.figure(7)
# plt.plot(freq_values, doppler_fft)
# plt.show()

# # dechirped_symbol = np.abs(np.fft.fft(symbol_train*doppler_out*DC, npt*len(DC[0:int(len(symbol_train)/2)])))
# # plt.figure(1)
# # plt.plot(dechirped_symbol[np.arange(-20,20)])
# # plt.show()
# # fp = FFO_corr_sym(symbol_train*doppler_out, SF, BW, FS, N, 0)
# # print(fp)