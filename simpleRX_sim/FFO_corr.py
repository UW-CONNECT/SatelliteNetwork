import numpy as np
from util import length
from stft_v2 import stft_v2
# from param_configs import param_configs
# from css_demod_2 import sym_to_css
import math, cmath
from matplotlib import pyplot as plt
import scipy.io
import scipy.signal

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

def FFO_corr(Data_stack,Upchirp_ind, SF, BW, Fs, N, num_preamble):
    '''
    Finds sampling offset for an upsampled signal Fs/BW that most 
    accurately dechirps the signal (i.e., minimizes leakage into adjacent bins
    and maximizes the peak)
    '''
    upsampling_factor = int(Fs/BW)
    
    DC = np.conjugate(sym_2_css([0], N, SF, BW, Fs))
    
    inn = []
    for i in range(0,10):
        freq_mean = [] 
        for pp in range(num_preamble-1):
            data_wind = Data_stack[int(Upchirp_ind+i)+(pp*N*upsampling_factor): int((Upchirp_ind+i)+((pp+1)*N*upsampling_factor)  )] 
            dech = np.abs(np.fft.fft(np.concatenate([(data_wind * DC), np.zeros((len(data_wind)))]), 32*len(data_wind)))
            dech = np.concatenate([dech[:int(N/2)], dech[-int(N/2):]])
            freq_track_qual = np.argmax(np.abs(dech))
            if (freq_track_qual > int(N/2)):
                freq_track_qual = N - freq_track_qual
            freq_mean.append(freq_track_qual**2) 
        freq_track_qual = np.mean(freq_mean)
        inn.append(freq_track_qual)
        
    inn = np.array(inn)
    b = inn.argmin(0) 
    return b 
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
    npt = 16
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
    # print("Tentative frequency offset: ", freq_offset)
        
    return -freq_offset 
    
    
def FFO_corr_sym_conj(Data_stack, SF, BW, Fs, N, sym):
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
    # true_sym = np.conjugate(sym_2_css([sym], N, SF, BW, Fs))
    true_sym = sym_2_css([sym], N, SF, BW, Fs)
    # npt = 32
    npt = 16
    
    # [TODO] Solve for base 2 number that provides the appropriate frequency resolution
    
    inn = []
    bin_vals = []
    freq_offsets = [] 
    # for i in np.arange(-5,5):
    for i in [0]:
        # freq_mean = [] 
        # for pp in range(num_preamble-1):
        data_wind = np.roll(Data_stack, i) # simulate different starting frequencies 
        # data_wind = np.conjugate(data_wind)
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
        bin_vals.append(-bin_val) 
        if (bin_offset > npt):
            # print("Positive offset.")
            bin_offset =  bin_offset - npt 
        else:
            bin_offset = -(npt - bin_offset)
        freq_offset = (Fs) * bin_offset/(N*upsampling_factor*npt) #* 10
        freq_offsets.append(freq_offset) 
        # print(bin_offset, "OFFSET")
    freq_offset = freq_offsets[np.argmax(bin_vals)]
    # print("Tentative frequency offset: ", freq_offset)
        
    return -freq_offset 