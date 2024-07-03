# CSS for Ham Satellites
Experimental 2500 Hz Bandwidth CSS for Ham Radio Linear Satellites (UHF/VHF) 
Ham operator call sign: KJ5DKX

# Instructions for Use 
config.py -> Modify these parameters to match your preferences (file location, chirp characteristics, etc) 

simpleTX_sim/genExperimentFile -> Generates CSS chirp samples and saves them to a file along with 
ground truth information. 

simpleRX_plusDoppler_simulation -> Change the file source to exp_root_folder/exp_folder/trial1, 
adjust SNR by changing either the multiply const or Channel Model 
[This is for simulation. Otherwise, just transmit trial1 over-the-air using USRP sink directly from TxfromFile.grc
and receive with simpleRx.grc]
 
simpleRX_sim/rx_stream_experiment.py -> Run, wait for "Started" message, then start the GNURADIO simulation. 

If the file paths and parameters in config.py were set correctly, packet information should begin to appear 
in the cmd window. Packet statistics such as PER and SER and doppler can be evaluated via [] 