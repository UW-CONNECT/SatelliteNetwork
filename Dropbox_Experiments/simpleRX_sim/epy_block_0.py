"""
Embedded Python Blocks:

Each time this file is saved, GRC will instantiate the first class it finds
to get ports and parameters of your block. The arguments to __init__  will
be the parameters. All of them are required to have default values!
"""

import numpy as np
from gnuradio import gr
import time 
import numpy as np 
import math
class blk(gr.sync_block): # other base classes are basic_block, decim_block, interp_block
    """Embedded Python Block example - a simple multiply const"""

    def __init__(self, sampling_rate=200e3, starting_frequency=0.0, frequency_slope=0.0):  # only default arguments here
        """arguments to this function show up as parameters in GRC"""
        gr.sync_block.__init__(
            self,
            name='Embedded Python Block',   # will show up in GRC
            in_sig=[np.complex64],
            out_sig=[np.complex64]
        )
        # if an attribute with the same name as a parameter is found,
        # a callback is registered (properties work, too).
        # self.example_param = example_param
        self.current_frequency = starting_frequency
        self.frequency_slope = frequency_slope
        self.start_time = time.time()
        self.sampling_rate = sampling_rate
        self.nSamp = 0 
        
    def work(self, input_items, output_items):
        """example: multiply with constant"""
        # current_time = time.time() - self.start_time 
        # print(current_time)
        current_time = self.nSamp / self.sampling_rate 
        # self.nSamp = self.nSamp + 1
        self.nSamp = self.nitems_written(0) / self.sampling_rate
        current_doppler = np.exp(1j *(2 * math.pi * self.frequency_slope*(current_time**2)/2))
        
        output_items[0][:] = input_items[0] * current_doppler
        return len(output_items[0])
