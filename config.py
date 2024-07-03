'''
config.py :: Jackie Schellberg :: Last updated 6/12/2024 


'''

''' General Params '''
RAW_FS = 200000     # Number of samples rx_stream waits from ZMQ until sending to demodulation 
FS = 200000         # (Hz) USRP Sampling Rate (FS/2 corresponding to the capture bandwidth) 
BW = 2500           # (Hz) 


''' Experiment Related ''' 
exp_root_folder = '../../experiment_data_JUNE2'                 # Toplevel file containing multiple experiments 
exp_folder = 'SF_7N_128BW_2500FS_200000NPKTS_50PLEN_100CR_0'    # Specific trial folder containing raw data and experiment details.
