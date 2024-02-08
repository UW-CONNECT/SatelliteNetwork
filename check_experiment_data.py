import pickle 
import glob 
import numpy as np

exp_root_folder = 'ber_desktop_testing'
exp_folder = '10pkts_sf9_20kbw_payload100'
trial_name = "trial1"
experiment_base_folder  = 'simpleTX_sim/'+exp_root_folder + '/' + exp_folder + '/' + 'error_out' + '*' + '.pkl'

files_list = glob.glob(experiment_base_folder)
print(files_list)
for file in files_list:
    with open(file,'rb') as f: 
        print(file)
        GND_TRUTH_PKT, OUTPUT_PKT, TOTAL_PKT_CNT = pickle.load(f)
        print(sum(np.subtract(GND_TRUTH_PKT,OUTPUT_PKT)))
        
