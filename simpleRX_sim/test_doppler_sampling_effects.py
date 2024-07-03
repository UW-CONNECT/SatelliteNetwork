''' Generate pure simulation chirps and measure the effects of frequency offset, sampling rate, etc '''
import math, cmath
from matplotlib import pyplot as plt
import scipy.io
import scipy.signal
import numpy as np
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

if __name__ == "__main__":
    symbol = [0]
    SF = 7 
    N = 2**SF
    BW = 2500
    FS = 200000
    # FS = 100000
    # FS = 50000
    doppler = 5003.5
    npt = 4
    tile_f  =2
    UPSAMP = int(FS/BW)
    sym_chirp = sym_2_css(symbol, N, SF, BW, FS) 
    
    tt = np.arange(0, len(sym_chirp))/FS
    sym_chirp = sym_chirp  * np.exp(1j * 2 * math.pi * doppler * tt)
    sym_chirp = np.tile(sym_chirp,tile_f)
    
    dc = np.conjugate(sym_2_css([0], N, SF, BW, FS))
    dc = np.tile(dc, tile_f)
    temp_wind_fft = abs(
            np.fft.fft(sym_chirp * dc, n=npt*len(dc), axis=0))        
    max_bin = np.argmax(temp_wind_fft)
    
    if (max_bin > ((N*UPSAMP*npt*tile_f))):    
        max_bin = -((N*UPSAMP*npt*tile_f) - max_bin)
    freq_shift = FS * max_bin/(N*UPSAMP*npt*tile_f)
    # print(freq_shift, "freq shift preabmle", "Freq resolution: ", self.FS * 1/(self.N*self.UPSAMP*npt), "samp len:", len(win1))
    possible_dop=int(freq_shift)
    
    print("Ground truth doppler: ", doppler, "Actual: ", possible_dop )
    print("Chirp time: ", len(sym_chirp)/FS)
    print("Waveform resolution: ", 1/(len(sym_chirp)/FS))
    print("Freq resolution: ", FS * 1/(N*UPSAMP*npt), "samp len:", len(sym_chirp))