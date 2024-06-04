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
                    
            #t_arr = np.arange(0, spsym*st_offset/N - 1)/Fs # [TODO] fix this ... 
            t_arr = np.arange(0, spsym*st_offset/N )/Fs 
            ce = np.exp(1j*(ang_start + 2*pi*(t_arr*(f_start +0.5*k*t_arr)) ))
            out_sig = np.concatenate((cs, ce))
            #sym_out.append(out_sig)
            sym_out=np.append(sym_out,out_sig)
           
        return sym_out

def FFO_corr(Data_stack,Upchirp_ind, SF, BW, Fs, N, num_preamble):
    upsampling_factor = int(Fs/BW)
    
    DC = np.conjugate(sym_2_css([0], N, SF, BW, Fs))
    
    inn = []
    # data_shifts = []
    # for i in range(upsampling_factor):
    for i in range(0,10):
    # for i in range(4):
        # freq_off = []
        # data_wind = Data_stack[int(Upchirp_ind+i): int((Upchirp_ind+i)+(num_preamble)*N*upsampling_factor )]   
        # Spec = stft_v2(data_wind, N*upsampling_factor, DC,0,0)
        # temp = []
        # freq_track_qual_before = ( np.sum((np.abs(Spec[0,:]) - np.abs(Spec[-1,:]))**2 ) + ( np.sum((np.abs(Spec[0,:]) - np.abs(Spec[1,:]))**2) ))**2
            # freq_track_qual_tmp.append(( np.sum(np.abs(Spec[0,:])) - np.sum(np.abs(Spec[-1,:])) )+np.sum(np.abs(Spec[0,:])) + ( np.sum(np.abs(Spec[0,:])) - np.sum(np.abs(Spec[1,:])) ))
        # freq_track_qual = np.mean(freq_track_qual_tmp)
        # freq_track_qual =  np.sum(np.abs(Spec[0,:])) 
        
        
        freq_mean = [] 
        for pp in range(num_preamble-1):
            data_wind = Data_stack[int(Upchirp_ind+i)+(pp*N*upsampling_factor): int((Upchirp_ind+i)+((pp+1)*N*upsampling_factor)  )] 
            # data_wind = data_wind * np.exp(1j * (- np.angle(data_wind[0])))   #tentatively adding leftover phase for random errors 
            # print("Phase angle correction:", - np.angle(data_wind[0]))
            # dech = np.abs(np.fft.fft(data_wind * DC, 32*len(data_wind)))
            dech = np.abs(np.fft.fft(np.concatenate([(data_wind * DC), np.zeros((len(data_wind)))]), 32*len(data_wind)))
            # fft_idx = np.concatenate([ np.arange(int(N/2)), np.arange(-int(N/2))
            dech = np.concatenate([dech[:int(N/2)], dech[-int(N/2):]])
            # freq_track_qual = ( np.sum((np.abs(dech[0]) - np.abs(dech[-1]))**2 ) + ( np.sum((np.abs(dech[0]) - np.abs(dech[1]))**2) ))**2
            freq_track_qual = np.argmax(np.abs(dech))
            if (freq_track_qual > int(N/2)):
                freq_track_qual = N - freq_track_qual
            # freq_track_qual = np.abs(dech[0])
            
            # Spec = stft_v2(data_wind, N*upsampling_factor, DC,0,0)
            # freq_track_qual = ( np.sum((np.abs(Spec[0,:]) - np.abs(Spec[-1,:]))**2 ) + ( np.sum((np.abs(Spec[0,:]) - np.abs(Spec[1,:]))**2) ))
            freq_mean.append(freq_track_qual**2) 
        freq_track_qual = np.mean(freq_mean)
        inn.append(freq_track_qual)
        
        # print("Freq track qual before: ", freq_track_qual_before, "After:", freq_track_qual)
    # plt.figure(1)
    # plt.plot(inn)
    # plt.show()
    inn = np.array(inn)
    # peak_stats.append(k_peak_stats)
    # Data_freq_off = np.array(Data_freq_off)
    # choosing the best Data_stack based on maximum energy difference from
    # adjacent bins
    # b = inn.argmax(0)
    b = inn.argmin(0) 
    # output frequency offset corrected buffer with relevant, Peak
    # statistics and frequency offsets
    # Data_buff.append(Data_freq_off[b,:])
    
    # Data_buff = Data_buff[b:]   # apply the sampling offset a sa frequency correction
    # FFO.append(ffo[b])
    # peak_amp.append(peak_stats[k][b])
    # Up_ind.append(Upchirp_ind[k,:])
    return b 