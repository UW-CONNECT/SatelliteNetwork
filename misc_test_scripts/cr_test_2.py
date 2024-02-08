import numpy as np 

def hamming_encode(nibbles):
    n_count = len(nibbles)
    codewords = np.zeros(n_count)
    
def bit_reduce(w,pos):
    '''
    b = bool((w&(2**pos[0])) > 0)
    for i in range(1,len(pos)):
        b = b ^ bool( (w & (2**pos[i])) > 0 )
        print(b,  w & (2**pos[i]), pos[i])
    '''
    xor_val = 0
    for i in range(0, len(pos)):
        xor_val = xor_val + 2**pos[i] 
    
    b = (((xor_val ^ w) ) == 0) 
    print(xor_val, (xor_val ^ w) )
    return b
        
if __name__ == "__main__":
    print(bit_reduce(14, [1,2,3,4]))