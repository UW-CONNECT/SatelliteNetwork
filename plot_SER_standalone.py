'''
For processing the experiment data into plots and such, given ground truth data. 

Requires a POSIX doppler start time, which should be copied and placed into the appropriate folder for the 
trial under test.
'''
import pickle 
import glob 
import numpy as np
from matplotlib import pyplot as plt
import math

 # = r'J:\schellberg\indoor_exp_feb_2024\SatelliteNetwork-main\simpleTX_sim\doppler_sim_scripts\doppler_sync_testing_time.pkl'
# # output_time_vec = np.arange(0, len(output))/ FS 
# f = opdoppler_time_fileen(doppler_time_file,'rb')
# [doppler_ref, FS_dop] = pickle.load(f)
# doppler_ref = np.array(doppler_ref)
# doppler_ref = doppler_ref[:,0]
# # doppler_ref = np.concatenate([doppler_ref, np.flip(doppler_ref)])
# f.close()

# FS_dop = 100
# doppler_time_vec = np.arange(0, len(doppler_ref))/FS_dop
# print(doppler_time_vec[-1],len(doppler_ref),len(doppler_ref)/FS_dop)
# plt.figure(1)
# plt.plot(doppler_time_vec, doppler_ref)
# plt.xlabel('Time (s)')
# plt.ylabel('Doppler Shift (Hz)')
# plt.title('Doppler Shift Simulated with Middle USRP (taken from an ISS pass)')
# plt.show()

# plt.figure(1)
# plt.hist(np.abs(np.diff(doppler_ref[::int(FS_dop/20)])))
# plt.show()
# print("Max doppler change: ", np.abs(np.diff(doppler_ref[::int(FS_dop/20)])))
SF = [7, 8, 9, 10, 11, 12]
BW = 20000
FS = 200000
NUMPKTS = 10 
PAYLOAD_LEN = 100
#exp_root_folder = 'indoor_testing'
# exp_root_folder = '../experiment_data_maxtx1transponder'
#exp_root_folder = '../simulated_data'
#exp_root_folder = '../simulated_data'
# exp_root_folder = '../experiment_data'
SF_x = []
SF_pkt_x = []
PKT_CNT_Y = []
SER_Y = []
expected_pkt_cnt = 10 # this should be saved in our ground truth file... 
#for sf in SF: 

# exp_root_folder = '../experiment_data_synced'
# exp_root_folder = '../experiment_data_june2'

exp_root_folder = '../EXPERIMENT_DATA_7_8_2024'

f = open('trial_under_test.txt')
test_list = [f.read()]
f.close()

# doppler Start time 
# f = open('J:\schellberg\indoor_exp_feb_2024\SatelliteNetwork-main\doppler_start_time_jun23')
# d_start_time = float(f.read())
# f.close()

# print("Doppler start time: ", d_start_time)

output_file_cnt = 0
SNR_dx = []
error_locs = []
rx_packet_time_vec = []
for exp_folder in test_list:
    sf=7
    experiment_base_folder = exp_root_folder + '/' + exp_folder + '/' + 'error_out' + '*' + '.pkl'
    files_list = glob.glob(experiment_base_folder)
    # print(len(files_list))
    
    output_file_cnt = 0
    for file in files_list:
        with open(file,'rb') as f: 
            #print(file)
            # GND_TRUTH_PKT, OUTPUT_PKT, TOTAL_PKT_CNT, PKT_SNR = pickle.load(f)
            GND_TRUTH_PKT, OUTPUT_PKT, TOTAL_PKT_CNT, PKT_SNR,CURR_POSIX_TIME = pickle.load(f)
            # rx_packet_time_vec.append(CURR_POSIX_TIME)
            # CURR_POSIX_TIME = CURR_POSIX_TIME % d_start_time 
            rx_packet_time_vec.append(CURR_POSIX_TIME)
            if len(OUTPUT_PKT) == len(GND_TRUTH_PKT): 
                SF_x.append(sf)                
                n_errors = int(np.count_nonzero(abs(np.subtract(GND_TRUTH_PKT,OUTPUT_PKT))) ) 
                if (n_errors > 0):
                    # print("DX",np.argwhere(abs(np.subtract(GND_TRUTH_PKT,OUTPUT_PKT)) > 0))
                    error_locs.extend(np.argwhere(abs(np.subtract(GND_TRUTH_PKT,OUTPUT_PKT)) > 0))
                    # plt.figure(2)
                    # plt.plot(abs(np.subtract(GND_TRUTH_PKT,OUTPUT_PKT)))
                    # plt.show()
                SER_Y.append(n_errors) 
                output_file_cnt += 1
                SNR_dx.append(PKT_SNR)
            else:
                SER_Y.append(100)
    SF_pkt_x.append(sf)
    PKT_CNT_Y.append(output_file_cnt) 


# find the corresponding ground truth doppler times, presuming it is correct. 
sorted_rx_time = np.argsort(rx_packet_time_vec) 
rx_packet_time_vec  = np.array(rx_packet_time_vec)
SER_Y = np.array(SER_Y)
rx_packet_time_vec = rx_packet_time_vec[sorted_rx_time]
SER_Y = SER_Y[sorted_rx_time]
# rx_packet_time_vec_old = rx_packet_time_vec
# rx_packet_time_vec = np.subtract(rx_packet_time_vec , rx_packet_time_vec[0])
# plt.figure(1)
# plt.plot(rx_packet_time_vec)
# plt.show()
# print("Total packets received:", output_file_cnt)
print("Total packets received:", output_file_cnt)
print("Average SNR:", np.mean(SNR_dx))
print("Mean SER", np.mean(SER_Y/len(GND_TRUTH_PKT)))
# # doppler_indices = []
# # for pt in rx_packet_time_vec: 
    # # doppler_indices.append(np.argmin(abs( np.subtract(doppler_time_vec,pt))))
# # dop_diff = abs(np.diff(doppler_ref)) 
# x = dop_diff[doppler_indices] 
# # x = doppler_ref[doppler_indices] 

plt.figure(1)
plt.scatter(SNR_dx, SER_Y)
plt.xlabel('SNR (dB)')
plt.ylabel('Number of Incorrect Symbols')

# plt.figure(2)
# plt.plot(rx_packet_time_vec,x)
# plt.xlabel('Time (s)')
# plt.ylabel('Doppler Shift (Hz)')
plt.show()

nandx = np.argwhere(np.isnan(SNR_dx))
