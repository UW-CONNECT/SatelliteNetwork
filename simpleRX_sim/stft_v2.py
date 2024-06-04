import numpy as np
import math
from matplotlib import pyplot as plt
def stft_v2(Rx_Buffer,N,DC,upsamp,dis):
    # This function produces a spectrogram of LoRa signal using Dechirping
    # operation to get the best frequency Resolution
    Spec = np.zeros((N,len(Rx_Buffer)))
    # Spec = np.zeros((128,len(Rx_Buffer)))
    buff = np.concatenate([Rx_Buffer, np.zeros(N-1)])
    if(upsamp):
        for i in range(len(Rx_Buffer)):
            # Spec[:,i] = np.roll(np.abs(np.fft.fft(buff[i:i+N] * DC.conj())) / math.sqrt(N),-round( (i)/8 ))
            Spec[:,i] = np.roll(np.abs(np.fft.fft(buff[i:i+N] * DC)) / math.sqrt(N),round( (i)/10 ))
    else:
        for i in range(len(Rx_Buffer)):            
            # Spec[:,i] = np.roll(np.abs(np.fft.fft(buff_tmp * DC)) / math.sqrt(N), (i))
            
            Spec[:,i] = np.roll(np.abs(np.fft.fft(buff[i:i+N] * DC)) / math.sqrt(N), (i))
            
            # Spec[:,i] = np.roll(np.abs(np.fft.fft(buff[i:i+N] * DC, 128)) / math.sqrt(N), (i))
            # plt.figure(1)
            # plt.plot(np.roll(np.abs(np.fft.fft(buff[i:i+N] * DC)) / math.sqrt(N), (i)))
            # plt.show()         

    return Spec

