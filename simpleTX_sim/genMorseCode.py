from css_mod import CssMod
from matplotlib import pyplot as plt
import numpy as np
import time
import math
import zmq
import random
from pathlib import Path
import pickle 

# css_modulator = CssMod(N, SF, BW, FS, preamble, end_delimeter, CR)
def string2cw(string):
    out_data = []
    for character in string.upper():
        morse_sym = morse_dict[character] 
        for morse in morse_sym: 
            if morse == '.':
                out_data.extend(short_sym)
            elif morse == '-':
                out_data.extend(long_sym)
            else: # new word 
                out_data.extend(word_space_sym)
            out_data.extend(intra_space_sym)
        out_data.extend(letter_space_sym)
    return out_data 

def ang2bin_nopad(data):
    '''
    Convert complex data to binary for GNUradio 
    '''
    Rx_buff = np.zeros((2*len(data)))
    Rx_buff[0::2] = np.real(data)
    Rx_buff[1::2] = np.imag(data)
   
    #Rx_buff = np.append(np.zeros((len(Rx_buff))),Rx_buff)
    
    return Rx_buff

if __name__ == '__main__':    
    FS = 200000 # typically, dot is 1 second and dash is 3x dot, space is 3x dot
    
    SHORT = FS
    LONG = SHORT * 3 
    
    INTRA_LETTER_SPACE = SHORT 
    LETTER_SPACE = FS * 3 
    WORD_SPACE = FS * 7 
    
    long_sym = ang2bin_nopad(np.ones(LONG))
    short_sym = ang2bin_nopad(np.ones(SHORT))
    intra_space_sym = ang2bin_nopad(np.zeros(INTRA_LETTER_SPACE))
    letter_space_sym = ang2bin_nopad(np.zeros(LETTER_SPACE))
    word_space_sym = ang2bin_nopad(np.zeros(WORD_SPACE))

    print(len(short_sym),len(intra_space_sym),len(letter_space_sym),len(word_space_sym))
    
    # dict from https://www.educative.io/answers/how-to-write-a-morse-code-translator-in-python 
    morse_dict = {
        'A': '.-', 'B': '-...', 'C': '-.-.', 'D': '-..', 'E': '.', 'F': '..-.', 'G': '--.', 'H': '....',
        'I': '..', 'J': '.---', 'K': '-.-', 'L': '.-..', 'M': '--', 'N': '-.', 'O': '---', 'P': '.--.',
        'Q': '--.-', 'R': '.-.', 'S': '...', 'T': '-', 'U': '..-', 'V': '...-', 'W': '.--', 'X': '-..-',
        'Y': '-.--', 'Z': '--..', '0': '-----', '1': '.----', '2': '..---', '3': '...--', '4': '....-',
        '5': '.....', '6': '-....', '7': '--...', '8': '---..', '9': '----.', '.': '.-.-.-', ',': '--..--',
        '?': '..--..', "'": '.----.', '!': '-.-.--', '/': '-..-.', '(': '-.--.', ')': '-.--.-', '&': '.-...',
        ':': '---...', ';': '-.-.-.', '=': '-...-', '+': '.-.-.', '-': '-....-', '_': '..--.-', '"': '.-..-.',
        '$': '...-..-', '@': '.--.-.', ' ': '/'
    }

    TX_STRING = 'KJ5DKX EN53hb' # Callsign and gridlocator for a basic contact 
    cw = string2cw(TX_STRING)
    print(len(cw))
    bin_dat = np.float32(cw).tobytes()
    
    file = open(TX_STRING, 'bw')
    file.write(bin_dat)
    file.close()