import numpy as np
# from util import length
# from stft_v2 import stft_v2
# from param_configs import param_configs
# from css_demod_2 import sym_to_css
import math, cmath
from matplotlib import pyplot as plt
import scipy.io
import scipy.signal
from stft_v2 import stft_v2
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
def FFO_corr_UPSAMP(Data_stack,Upchirp_ind, SF, BW, Fs, N, num_preamble):
    '''
    Finds sampling offset for an upsampled signal Fs/BW that most 
    accurately dechirps the signal (i.e., minimizes leakage into adjacent bins
    and maximizes the peak)
    '''
    upsampling_factor = int(Fs/BW)
    
    DC = np.conjugate(sym_2_css([0], N, SF, BW, Fs))
    
    inn = []
    # shifts = range(-N*upsampling_factor,N*upsampling_factor)
    # shifts = range(-upsampling_factor,upsampling_factor)
    shifts = np.arange(-upsampling_factor,upsampling_factor)
    # shifts = np.arange(-int(upsampling_factor/2),int(upsampling_factor/2))
    # shifts = range(0,upsampling_factor)
    for i in shifts:
        freq_mean = [] 
        for pp in [0]:
            data_wind = Data_stack[int(Upchirp_ind+i)+(pp*N*upsampling_factor): int((Upchirp_ind+i)+((pp+1)*N*upsampling_factor)  )] 
            # dech = np.abs(np.fft.fft(np.concatenate([(data_wind * DC), np.zeros((len(data_wind)))]), 32*len(data_wind)))
            # dech = np.abs(np.fft.fft(np.concatenate([(data_wind * DC), np.zeros((len(data_wind)))])))
            # # dech = np.abs(np.fft.fft((data_wind * DC)))
            # # # dech = np.concatenate([dech[:int(N/2)], dech[-int(N/2):]])
            # # # freq_track_qual = np.argmax(np.abs(dech))
            # # freq_track_qual_idx = np.argmax(np.abs(dech))
            # # # if (freq_track_qual_idx > int(N/2)):
                # # # freq_track_qual_idx = N - freq_track_qual_idx
            # # # freq_mean.append(freq_track_qual**2) 
            
            # # if (freq_track_qual_idx == len(dech)-1):
                # # s_dx = -1 
            # # else:
                # # s_dx = freq_track_qual_idx+1
            
            # # freq_mean.append((dech[freq_track_qual_idx]**2 - (dech[s_dx]+dech[freq_track_qual_idx-1])**2)**2) 
            # # print(freq_track_qual_idx, s_dx,freq_track_qual_idx-1, 'spec')
            print("Trying FFT")
            Spec = stft_v2(data_wind,N*upsampling_factor,DC,0,0)
            print("Done finding STFT")
            temp = []
            freq_track_qual = []
            pream_peak_ind = []
            adj_ind = []
            # row_ind contains the Frequency Rows around bin 1 where a
            # Preamble Peak can lie
            row_ind = np.concatenate([range(N-6,N), range(0,6)])
            count = 1
            for i in np.nditer(row_ind):
                temp.append(np.sum(np.abs(Spec[i,:])))
                count = count + 1
            temp = np.array(temp)
            # Frequency Track in row containing Preamble should have
            # maximum energy
            ind = temp.argmax(0)
            pream_peak_ind = row_ind[ind]
            # Find row indices for Preamble row + 1 & - 1
            adj_ind = np.array([np.mod(pream_peak_ind-1+1,N), np.mod(pream_peak_ind+1+1,N)]) # plus 1 for index conversion
            if(np.sum(adj_ind == 0) == 1):
                adj_ind[(adj_ind == 0).nonzero()] = N
            # A good quality frequency track for a preamble is one that has
            # least energy leakage in adjacent rows (this promises very sharp FFT peaks)
            adj_ind -= 1 # subtract 1 to convert back to Python indices
            freq_track_qual = ( np.sum(np.abs(Spec[pream_peak_ind,:])) - np.sum(np.abs(Spec[adj_ind[0],:])) ) + ( np.sum(np.abs(Spec[pream_peak_ind,:])) - np.sum(np.abs(Spec[adj_ind[1],:])) )
            freq_mean.append(freq_track_qual)
        freq_track_qual = np.mean(freq_mean)
        inn.append(freq_track_qual)
        
    inn = np.array(inn)
        
    plt.figure(7)
    plt.plot(shifts,inn)
    plt.title('')
    plt.show()
    
    # b = inn.argmin(0)
    # b = shifts[inn.argmin(0)] + 1
    b = shifts[inn.argmax(0)] #+ 1
    return b 
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
# def FFO_corr_sym(Data_stack, SF, BW, Fs, N, sym):
    # '''
    # Determines the fine frequency offset from the received samples in Datastack to the 
    # ideal generated symbol by dechirping with an upsampled FFT with resolution Fs/M and
    # identifying a frequency offset. 
    # '''
    # upsampling_factor = int(Fs/BW)
    # # print(Fs)
    # # print(upsampling_factor)
    # # DC = np.conjugate(sym_2_css([0], N, SF, BW, Fs))
    # # true_sym = sym_2_css([sym], N, SF, BW, Fs)
    # true_sym = sym_2_css([sym], N, SF, BW, Fs)
    # # npt = 16
    # # npt = 4#*2*2#*16
    # npt =1 
    # # npt = 2
    # # npt = 4
    
    # # [TODO] Solve for base 2 number that provides the appropriate frequency resolution
    # ts = np.arange(0, len(true_sym)) / Fs 
    # inn = []
    # bin_vals = []
    # freq_offsets = [] 
    # # for i in np.arange(-5,5):
    # possible_freqs = np.arange(-30,30,1) #[TODO] Replace this.     
    
    # corr_arr = []
    # for freq1 in possible_freqs: 
        # freq = np.linspace(0, freq1, len(true_sym) )   
        # corr = np.abs(np.correlate(Data_stack, true_sym *  np.exp(1j * 2 * math.pi *freq * ts)))
        # corr_arr.append(corr) 
    
    # plt.figure(1)
    # plt.plot(possible_freqs, corr_arr) 
    # plt.show() 
    
    # freq_offset = possible_freqs[np.argmax(corr_arr)]
       
    # return -freq_offset   
def FFO_corr_sym(Data_stack, SF, BW, Fs, N, sym):
    '''
    Determines the fine frequency offset from the received samples in Datastack to the 
    ideal generated symbol by dechirping with an upsampled FFT with resolution Fs/M and
    identifying a frequency offset. 
    '''
    upsampling_factor = int((Fs/BW )) 
    # print(Fs)
    # print(upsampling_factor)
    # DC = np.conjugate(sym_2_css([0], N, SF, BW, Fs))
    # true_sym = sym_2_css([sym], N, SF, BW, Fs)
    true_sym = np.conjugate(sym_2_css([sym], N, SF, BW, Fs))
    
    # Data_stack = Data_stack[::upsampling_factor]
    # true_sym = true_sym[::upsampling_factor]
    
    # npt = 16
    # npt = 4#*2*2#*16
    # npt =8 
    npt = 1
    # npt = 4
    
    # [TODO] Solve for base 2 number that provides the appropriate frequency resolution
    ts = np.arange(0, len(true_sym)) / Fs
    # ts = np.arange(0, len(true_sym)) / (Fs/upsampling_factor)    
    inn = []
    bin_vals = []
    freq_offsets = [] 
    # for i in np.arange(-5,5):
    possible_freqs = np.arange(-6,6,.1) #[TODO] Replace this. 
    max_val = []
    st_freqs = []
    end_freqs = []
    possible_start_freqs = np.arange(-6,6,.1) 
    possible_end_freqs = np.arange(-6,6,.1)
    for st_freq in possible_start_freqs:
        for freq1 in possible_end_freqs:
            # freq = (np.arange(0, len(true_sym)) / Fs ) * freq1 # ending frequency matters? 
            # freq = np.arange(0, freq1, freq1/len(true_sym) )
            # freq = np.linspace(0, freq1, len(true_sym) )  
            # freq = np.linspace(0, freq1, len(true_sym) ) 
            freq = np.linspace(st_freq, freq1, len(true_sym) )             

            # plt.figure(1)
            # plt.plot(freq)
            # plt.title("Rd")
            # plt.show()             
                   
            doppler_temp  = np.exp(1j * 2 * math.pi * -freq * ts)
            # doppler_temp  = np.conjugate(np.exp(1j * 2 * math.pi * freq * ts))
            # doppler_temp  = np.conjugate(np.exp(1j * 2 * math.pi * freq * ts))
            # freq_mean = [] 
            # for pp in range(num_preamble-1):
            # data_wind = np.roll(Data_stack, i) # simulate different starting frequencies 
            # dech = np.abs(np.fft.fft(np.concatenate([(data_wind * true_sym))))
            # dech1 = np.abs(np.fft.fft(np.concatenate([(data_wind * true_sym)]), npt*len(data_wind)))
            dech1 = np.abs(np.fft.fft((Data_stack * true_sym * doppler_temp), npt*N*upsampling_factor))
            
            
            # plt.figure(60)
            # plt.plot(np.abs(np.fft.fft((Data_stack * true_sym * doppler_temp), npt*N*upsampling_factor))[np.arange(-10,10)])
            # plt.title("Using Line ")
            # plt.figure(59)
            # plt.plot(np.abs(np.fft.fft((Data_stack * true_sym * freq1), npt*N*upsampling_factor)) [np.arange(-10,10)])
            # plt.title("Using Constant tone")
            # print("Frequency:", freq1)
            # plt.show()
            st_freqs.append(st_freq)
            end_freqs.append(freq1)
            max_val.append(dech1[0]**2 - (dech1[-1]**2 +dech1[1]**2)) 
            # max_val.append(dech1[0]**2 - (dech1[-1] +dech1[1])**2) 
            
        # max_val.append(max(dech1)) 
    # possible_freqs = possible_freqs * len(true_sym) / Fs
    # freq_offset = possible_freqs[np.argmax(max_val)] 
    # print("Tentative frequency offset: ", freq_offset)
    
    # plt.figure(6)
    # plt.specgram(Data_stack*true_sym, Fs=Fs)
        # # plt.show()
    # plt.figure(60)
    # freq = np.linspace(0, freq_offset, len(true_sym) )      

    # # plt.figure(1)
    # # plt.plot(freq)
    # # plt.title("Rd")
    # # plt.show()
    
    # doppler_temp  = np.exp(1j * 2 * math.pi * -freq * ts)
    
    # plt.plot(np.abs(np.fft.fft((Data_stack * true_sym * doppler_temp), npt*N*upsampling_factor))[np.arange(-10,10)])
    # plt.title("Using Line ")
    # plt.figure(59)
    # doppler_temp  = np.exp(1j * 2 * math.pi * -freq_offset * ts)
    # plt.plot(np.abs(np.fft.fft((Data_stack * true_sym * doppler_temp), npt*N*upsampling_factor)) [np.arange(-10,10)])
    # plt.title("Using Constant tone")
    # print("Frequency:", freq_offset, "Time: ", ts[-1])
    # plt.show()
            
    # plt.figure(1)
    # plt.plot(st_freqs, max_val)
    # plt.show()
    # freq_offset = 10.1     
    
    max_dx = np.argmax(max_val) 
    slope = (end_freqs[max_dx] - st_freqs[max_dx]) #/ (ts[-1])  
    print("==Slope: (Hz/s)", slope/ (ts[-1]) , "Starting frequency :", st_freqs[max_dx], "Ending frequency: ", end_freqs[max_dx]) 
    # freq_offset = slope - st_freqs[max_dx]
    freq_offset = 2*slope #- st_freqs[max_dx]
    # freq_offset = slope + st_freqs[max_dx]
    
    print("Frequency:", freq_offset, "Time: ", ts[-1])   
    return -freq_offset     
    
# def FFO_corr_sym(Data_stack, SF, BW, Fs, N, sym):
    # '''
    # Determines the fine frequency offset from the received samples in Datastack to the 
    # ideal generated symbol by dechirping with an upsampled FFT with resolution Fs/M and
    # identifying a frequency offset. 
    # '''
    # upsampling_factor = int(Fs/BW)
    # # print(Fs)
    # # print(upsampling_factor)
    # # DC = np.conjugate(sym_2_css([0], N, SF, BW, Fs))
    # # true_sym = sym_2_css([sym], N, SF, BW, Fs)
    # true_sym = np.conjugate(sym_2_css([sym], N, SF, BW, Fs))
    # npt = 16
    # # npt = 64
    
    # # [TODO] Solve for base 2 number that provides the appropriate frequency resolution
    
    # inn = []
    # bin_vals = []
    # freq_offsets = [] 
    # # for i in np.arange(-5,5):
    # for i in [0]:
        # # freq_mean = [] 
        # # for pp in range(num_preamble-1):
        # data_wind = np.roll(Data_stack, i) # simulate different starting frequencies 
        # # dech = np.abs(np.fft.fft(np.concatenate([(data_wind * true_sym))))
        # # dech1 = np.abs(np.fft.fft(np.concatenate([(data_wind * true_sym)]), npt*len(data_wind)))
        # dech1 = np.abs(np.fft.fft(np.concatenate([(data_wind * true_sym)]), npt*N*upsampling_factor))
        # # dech1 = np.concatenate([dech1[:npt], dech1[-npt:]])
        # dech1 = np.concatenate([dech1[-npt:],dech1[:npt]])
        # # plt.figure(2)
        # # plt.plot(dech1)
        # # plt.show()
        
        # bin_offset = np.argmax(dech1)
        # bin_val = np.max(dech1) 
        # bin_vals.append(bin_val) 
        # if (bin_offset > npt):
            # # print("Positive offset.")
            # bin_offset =  bin_offset - npt 
        # else:
            # bin_offset = (npt - bin_offset)
        # freq_offset = (Fs) * bin_offset/(N*upsampling_factor*npt) #* 10
        # # print("Resolution: ", (Fs) * 1/(N*upsampling_factor*npt))
        # freq_offsets.append(freq_offset) 
        # # print(bin_offset, "OFFSET")
    # freq_offset = freq_offsets[np.argmax(bin_vals)]
    # # print("Tentative frequency offset: ", freq_offset)
        
    # return -freq_offset 
    
    
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
    # npt = 16
    npt = 4
    
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