"""
Embedded Python Blocks:

Each time this file is saved, GRC will instantiate the first class it finds
to get ports and parameters of your block. The arguments to __init__  will
be the parameters. All of them are required to have default values!
"""

import numpy as np
from gnuradio import gr
import time 
import struct
class blk(gr.sync_block):  # other base classes are basic_block, decim_block, interp_block
    """Embedded Python Block example - a simple multiply const"""

    def __init__(self,  doppler_file_name=r'J:\schellberg\indoor_exp_feb_2024\doppler_simulation_files\GNURADIO_linear_8k_sep'):  # only default arguments here
        """arguments to this function show up as parameters in GRC"""
        gr.sync_block.__init__(
            self,
            name='Apply Doppler and Save Posix Time',   # will show up in GRC
            in_sig=[],
            out_sig=[np.complex64]
        )
        # if an attribute with the same name as a parameter is found,
        # a callback is registered (properties work, too).
        
        # self.example_param = example_param
        self.doppler_file_name = doppler_file_name
        
    
        # experiment_time_path= EXP_ROOT_FOLDER + '\\' + str(trial_number) + '\\' + "posix_start_time"
        self.f = open(doppler_file_name, "rb")
        # self.doppler_x = np.frombuffer(self.f.read(), dtype=np.float32)
        # self.doppler = self.doppler_x[::2] + 1j*self.doppler_x[1::2]
        # # with open(doppler_file_name, "rb") as f:
            # # print("Opened.")
            # nbytes= int((32 /8) *2 * 200000)
            # # nbytes= 8
            # # while (byte := f.read(nbytes)):
                # Do stuff with byte.        
                # zsocket.send(byte)
                # byte = np.frombuffer(byte, dtype=np.float32)
                # print("Sent byte")
                # time.sleep(.1)\
                # # print(byte)
            # self.doppler = struct.unpack('f', f.read())
        # print(len(self.doppler), "DOpple")
        # f.close()
        self.doppler_ptr = 0
        experiment_start_time = time.time()
        print("Experiment Start Time:",experiment_start_time)
        
    def work(self, input_items, output_items):
        doppler_re = np.complex64(struct.unpack('f',self.f.read(4)))
        doppler_im = np.complex64(struct.unpack('f',self.f.read(4)))
        # [doppler_re, doppler_im]=
        # doppler_im=np.frombuffer(self.f.read(4), dtype=np.float32)
        # print(doppler_re, doppler_im)
        # doppler_im=np.frombuffer(self.f.read(4), dtype=np.float32)
        doppler = doppler_re + 1j*doppler_im 
        # print(doppler)
        # doppler = 1
        # if (self.doppler_ptr == len(self.doppler)):
            # self.doppler_ptr = 0 
        # else: 
            # self.doppler_ptr = self.doppler_ptr + 1
        # output_items[0][:] = input_items[0] * self.doppler[self.doppler_ptr]
        output_items[0][:] = doppler
        
        return len(output_items[0])
    def stop(self):
        self.f.close()
        return True