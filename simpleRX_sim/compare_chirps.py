from css_demod_2 import CssDemod
import numpy as np
SF = 9
N = 2**SF 
UPSAMP = 10
preamble = [1,1,1]
PREAMBLE_SZ = len(preamble) 
END_DELIMITER = [1,1,1,1] 
symbols = np.array([1,1,1]) 
PAYLOAD_LEN = 100 
NUMPKTS = 10 
DB_THRESH = -8
css_demodulator = CssDemod(N, UPSAMP,PREAMBLE_SZ,END_DELIMITER,DB_THRESH, symbols,PAYLOAD_LEN,NUMPKTS,SF);

#css_demodulator.scipy_chirp_shift([1],N, UPSAMP)
print("Starting chirp comparison")
#(self, symbol, N, SF, BW, Fs)
BW = 20000
Fs = 200000
css_demodulator.sym_2_css([1,50,400],N, SF,BW,Fs)
