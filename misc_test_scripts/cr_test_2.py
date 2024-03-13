'''
cr_test_2.py :: Jackie Schellberg :: Last Updated Feb 2024 

For testing hamming encode and decode. 


Based off of: 
https://github.com/fabio-busacca/sdr-lora 
https://github.com/jkadbear/LoRaPHY/blob/master/README.md
'''

import numpy as np 

def hamming_encode(nibbles, CR):
    '''
    These are setup to take 8-bit ints. May need to change this later. 
    
    nibbles = 4 bit 
    CR = coding rate (1-4) 
    '''
    n_count = len(nibbles)
    #codewords = np.zeros(n_count)
    codewords = []
    
    for i in range(0, n_count):
        nib = nibbles[i] 
        p1 = bit_reduce(nib, [0, 2, 3])
        p2 = bit_reduce(nib, [0, 1, 3])
        p3 = bit_reduce(nib, [0, 1, 2])
        p4 = bit_reduce(nib, [0, 1, 2, 3])
        p5 = bit_reduce(nib, [1, 2, 3])
        
        if (CR == 1):
            codewords.append(p4 * 16 | nib)
        elif (CR == 2):
            #codewords.append(p5*32 )
            codewords.append(p5*32 | p3*16 | nib)
        elif (CR == 3):
            codewords.append(p2*64 | p5*32 | p3*16 | nib)
        elif (CR == 4):
            codewords.append(p1*128 | p2*64 | p5*32 | p3*16 | nib)
        else: 
            print("Check CR ... [encode]")
        
    return codewords 
    
def hamming_decode(codewords, CR):
    '''
    Hamming decode, where CR is the coding rate. 
    '''
    CR = CR + 4
    nibs = []
    print("Coding rate + 4: ", CR)
    for i in range(0, len(codewords)):
        codeword = codewords[i]
        p1 = bit_reduce(codeword, [7,3,2,0])
        p2 = bit_reduce(codeword,[6,3,1,0])
        p3 = bit_reduce(codeword,[4,2,1,0])
        p4 = bit_reduce(codeword,[4,3,2,1,0])
        p5 = bit_reduce(codeword,[5,3,2,1])
        
        print("Current codeword:", codeword) 
        
        # implement this 
        if CR == 5 or CR == 6:
            # detect but don't correct 
            nibs.append(np.mod(codeword, 16)) 
        elif CR == 7 or CR == 8: 
            par = p2*4 + p3*2 + p5 
            # fix the parity and modify par  
            #print("Parity:", par)
            
            if par == 3: 
                par = 4 
            elif par == 5:
                par = 8
            elif par == 6:
                par = 1 
            elif par == 7:
                par = 2
            else:
                par = 0
            
            codeword = codeword ^ par 
            nibs.append(np.mod(codeword, 16)) 
            # pass
        else:
            print("WRONGO BUDDY!!!")
    return nibs 
    #p1 = bit_reduce()
#def bit_or(w1,w2):
    

def bit_reduce(w,pos):
    '''
    XOR two integers -> Essentially 'and' w with pos, and count the number of 1's mod 2 
    
    0-based indexing for pos[i]
    '''
    b = 0
    #pos = pos - 1
    for i in range(0, len(pos)):
        if ((w & 2**(pos[i])) > 0):
            b = b + 1
    return np.mod(b,2)==1
        
if __name__ == "__main__":
    #print(bit_reduce(14, [1,2,3,4]))
    #print(bit_reduce(14, [2,3,4]))
    int_example = 0
    CR = 4
    # encoded = hamming_encode([int_example & 15, int_example & 240], CR)
    encoded = hamming_encode([int_example], CR)
    
    # Test a single bit error. 
    print("Encoded before", encoded)
    encoded = encoded[0] | 3            # 1, 2, 4, or 8  
    print("Encoded after: ", encoded) 
    
    decoded = hamming_decode([encoded], CR)
    print("Decoded:" ,decoded)