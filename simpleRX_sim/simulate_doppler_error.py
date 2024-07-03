import math
import multiprocessing
from multiprocessing import Manager
import socket

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
FS =200000
# FS =5000
SLOPE_TIME = 5*60*FS # number of samples, minutes * seconds * samples/second 
doppler_slope = 100 # Rate of change of doppler, Hz/second

# starting_frequency = 0
# starting_frequency = 0

time_vec = (np.arange(0, SLOPE_TIME) / FS )
# doppler = (np.arange(0, SLOPE_TIME) / FS ) * doppler_slope
# doppler = (np.arange(0, SLOPE_TIME) / FS ) * doppler_slope

# output = np.exp(1j*( 2*math.pi*(doppler*time_vec) ))
SF =7 
N = 2**SF 
BW = 2500 
BW = 20000
# FS
symbol_size = len(sym_2_css([0],N,SF,BW,FS))
# symbol_train = np.concatenate((sym_2_css([0,0,0,0,0,0,0,0],N,SF,BW,FS),np.conjugate(sym_2_css([0,0],N,SF,BW,FS))))
# symbol_train = np.concatenate((sym_2_css([0],N,SF,BW,FS),np.conjugate(sym_2_css([0],N,SF,BW,FS))))
# symbol_train = np.concatenate((sym_2_css([0],N,SF,BW,FS),np.conjugate(sym_2_css([0],N,SF,BW,FS))))
# symbol_train = np.concatenate((sym_2_css([54,54,54,0,0],N,SF,BW,FS),np.conjugate(sym_2_css([0,0,60,5,3],N,SF,BW,FS))))
# conj_train = np.conjugate(symbol_train)
symbol_train = sym_2_css([0],N,SF,BW,FS)
# conj_train = np.conjugate(symbol_train)
# symbol_train = np.concatenate((sym_2_css([0],N,SF,BW,FS),np.conjugate(sym_2_css([0],N,SF,BW,FS))))

plt.figure(1)
plt.specgram(symbol_train)

plt.figure(2)
ddr = np.fft.ifft( np.fft.fft(symbol_train)+ np.fft.fft(symbol_train))
# ddr = symbol_train * symbol_train 
plt.specgram(ddr)

plt.figure(3)
plt.specgram(np.conjugate(ddr) * symbol_train)
plt.show()

# time_vec = (np.arange(0, len(symbol_train)) / FS )


# # doppler = (np.arange(0, SLOPE_TIME) / FS ) * doppler_slope
# doppler = ((np.arange(0, len(symbol_train)) / FS ) * doppler_slope) -8000
# doppler = np.exp(1j*( 2*math.pi*(doppler*time_vec) ))


# # doppler_100hz = (np.arange(0, len(symbol_train)) / FS ) * 100
# # doppler_100hz = np.exp(1j*( 2*math.pi*(doppler_100hz*time_vec) ))
# # conj_train = np.conjugate(symbol_train * doppler_100hz) 
# conj_train = np.conjugate(symbol_train ) 
 # #+ np.conjugate(symbol_train )   + np.conjugate(symbol_train *-doppler)  
# symbol_train = np.roll(symbol_train, 0)
# # plt.figure(8)
# # plt.specgram(conj_train)
# # plt.show()

# # conj_train = np.conjugate(symbol_train + symbol_train*np.exp(1j*( 2*math.pi*(doppler*time_vec) )))
# max_val = []
# # shifts=np.arange(-100, 100,int(FS/BW))
# shifts=np.arange(-20*int(FS/BW), 20*int(FS/BW),int(FS/BW))
# # shifts=np.arange(-10, 20,5)
# for i in shifts:
# # for i in [0]:
    # symbol_train_t = np.roll(symbol_train, i)

    # # plt.figure(1)
    # # plt.specgram(symbol_train * doppler) 
    # # plt.title(i)
    # # plt.show()
    # # plt.figure(2)
    # # plt.specgram(symbol_train * doppler) 
    # # plt.figure(3)
    # ins_idx = np.arange(-100,100)
    
    # doppler_fft = np.abs(np.fft.fft(symbol_train_t * doppler * conj_train) )
    # max_val.append(max(doppler_fft))
    # # max_val.append(max(doppler_fft[np.arange(-10,10)]))
    # # max_val.append(20*np.log10(doppler_fft[0]))
    # # plt.plot(doppler_fft[ins_idx])
    # # plt.title(i)
    # # plt.show()

# npt = .5
# max_bin = shifts[np.argmax(max_val)]
# print("Initial shift: ", shifts[np.argmax(max_val)])
# symbol_train_dop = np.roll(symbol_train, max_bin)
# doppler_fft = np.abs(np.fft.fft(symbol_train_dop * doppler * conj_train, int(len(symbol_train_dop)*npt)))
# max_bin = np.argmax(doppler_fft)

# plt.figure(4)
# plt.plot(shifts, max_val)
# plt.show()
# # plt.figure(6)
# # plt.plot(np.arange(-100,100),doppler_fft[np.arange(-100,100)])
# # plt.show()
# if (max_bin > ((N*int(FS/BW)*npt)/2)):    
    # max_bin = -((N*int(FS/BW)*npt) - max_bin)
# freq_shift = FS * max_bin/(N*int(FS/BW)*npt)
# print(freq_shift)
# st_dx = len(symbol_train)
# symbol_train = np.concatenate((symbol_train, symbol_train, symbol_train))
# shift_small = FFO_corr_UPSAMP(symbol_train,st_dx, SF, BW, FS, N, 1)
# print("Shifts small: ", shift_small)

# plt.figure(1)
# plt.plot(shifts, max_val)
# plt.show()