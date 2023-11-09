'''
test_css_demod.py :: Jackie Schellberg :: 10/20/2023 

For testing css_demod() without streaming UDP samples. 
Import saved ds_x_1 variable from matlab: save('testing_CSS_RX.mat','ds_x_1', '-v7')
'''
import pickle
from matplotlib import pyplot as plt
import numpy as np
from css_demod import CssDemod
import math

import scipy.io
# just transmitted a long string of ones; this should be line @ y = 3
ds_x_1 = scipy.io.loadmat('testing_CSS_noPhaseAcc_RX.mat')
queue = ds_x_1['ds_x_1']
queue = np.squeeze(queue)
#queue = queue[25000:]

# CSS Specifics to be known ahead of time
SF = 7 
N = 2**SF
UPSAMP = 10;
css_demodulator = CssDemod(N, UPSAMP);

# induce a fake doppler shift 
freq_shift = 9000
#freq_shift = 8000
#freq_shift = 6000
#freq_shift = 1000
#freq_shift = 0
t = np.linspace(0, len(queue)/200000, len(queue))

plt.figure(2)
xpts = range(0,len(queue))
plt.plot(xpts, queue)

queue = queue * np.exp(1j * 2 * math.pi * freq_shift * t)

output = []
css_demodulator.css_demod([], queue, output) 

plt.figure(0)
xpts = range(0,len(output))
plt.plot(xpts, output)

plt.figure(1)
xpts = range(0,len(queue))
plt.plot(xpts, queue)
plt.show()  
