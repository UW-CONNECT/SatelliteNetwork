'''
For testing hamming coding rates on symbols, similar to lora 

Some sources to check out:
https://medium.com/nerd-for-tech/intro-to-convolutional-coding-part-ii-d289c109ff7a
--> Suggests something about the Viterbi algorithm 

'''

import numpy as np

def hamming_CR(SF,CR, symbols):
    '''
    Encodes symbols with Hamming Code (like LoRa) according to the Coding Rate. 
    
    SF - spreading factor
    CR - coding rate (1-4)
    symbols - symbols (16 bit int) to be transmitted) 
    '''
    H = [np.array(([1], [1], [1], [1]), dtype=np.uint8),
        np.array(([1, 0], [1, 1], [1, 1], [0, 1]), dtype=np.uint8),
        np.array(([1, 0, 1], [1, 1, 1], [1, 1, 0], [0, 1, 1]), dtype=np.uint8),
        np.array(([1, 0, 1, 1], [1, 1, 1, 0], [1, 1, 0, 1], [0, 1, 1, 1]), dtype=np.uint8)]
    hamming_code = H[CR-1]
    
    # convert to binary 
    bin_sim = np.array(np.zeros(len(symbols)*8))
    for i in range(0,len(symbols)):
        bin_sim[i*8:i*8+8] = num2binary(symbols[i], 8)
    
    # add hamming codes 
    bin_sim_cr = np.array(np.zeros(len(symbols)*8 + 2*CR*len(symbols)))
    for i in range(0,int(len(bin_sim)/4)):
        '''
        print(np.concatenate( ( bin_sim[i*4:i*4 + 4], np.mod(bin_sim[i*4:i*4 + 4] @ hamming_code, 2) ) ).shape)
        print(bin_sim[i*4:i*4 + 4].shape)
        print( np.mod(bin_sim[i*4:i*4 + 4] @ hamming_code, 2).shape)
        '''
        bin_sim_cr[i*(4+CR):i*(4+CR)+(4+CR)] = \
            np.concatenate( ( bin_sim[i*4:i*4 + 4], np.mod(bin_sim[i*4:i*4 + 4] @ hamming_code, 2) ) )
    
    print(len(bin_sim_cr)/8)
    # covert back to symbols to modulate ==> note that bit lengths are SF bits long! 
    output_symbols = np.array(np.zeros(int(len(bin_sim_cr)/8)))
    for i in range(0, int(len(bin_sim_cr)/8)):
        output_symbols[i] = np.packbits( np.array(bin_sim_cr[i*8:i*8+8], dtype=bool))
    
    return output_symbols


def hamming_decode_CR():
    # hamming parity check matrices 
    
    # fix the parity error 
    total_bits = CR+4 
    if (total_bits == 5 || total_bits == 6):
        nibs = np.mod(codewords, 16) # check why this is the case 
    elif (total_bits == 7 || total_bits == 8):
        parity = p2*3 + p3*2 + p5
        
    else: 
        print("Wrong number of bits")
    pass



def binary2num(bits):
    num = np.packbits(bits)

def num2binary(num,length = 0):

    num = np.array([num],dtype=np.uint16)
    num = np.flip(num.view(np.uint8))
    num = np.unpackbits(num)
    return num[-length:]

if __name__ == "__main__":
    print(hamming_CR(SF=7,CR=1, symbols=[1]))