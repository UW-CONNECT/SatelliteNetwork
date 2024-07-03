'''
For processing the experiment data into plots and such, given ground truth data. 
'''
import pickle 
import glob 
import numpy as np
from matplotlib import pyplot as plt
import math

SF = [7, 8, 9, 10, 11, 12]
BW = 20000
FS = 200000
NUMPKTS = 10 
PAYLOAD_LEN = 100
#exp_root_folder = 'indoor_testing'
# exp_root_folder = '../experiment_data_maxtx1transponder'
#exp_root_folder = '../simulated_data'
#exp_root_folder = '../simulated_data'
exp_root_folder = '../experiment_data_JUNE17'
SF_x = []
SF_pkt_x = []
PKT_CNT_Y = []
SER_Y = []
expected_pkt_cnt = 10 # this should be saved in our ground truth file... 
#for sf in SF: 
# test_list = ['SF_7N_128BW_20000FS_200000NPKTS_10PLEN_100_0dBattn',
             # 'SF_7N_128BW_20000FS_200000NPKTS_10PLEN_100_20dBattn',
             # 'SF_7N_128BW_20000FS_200000NPKTS_10PLEN_100_50dBattn',
             # 'SF_7N_128BW_20000FS_200000NPKTS_10PLEN_100_2ndfloor']
# test_list = ['SF_7N_128BW_20000FS_200000NPKTS_10PLEN_100_0dBattn']
# test_list = ['SF_7N_128BW_20000FS_200000NPKTS_10PLEN_100_indoor7pwr']
# test_list = ['SF_7N_128BW_20000FS_200000NPKTS_10PLEN_100']
# test_list = ['SF_7N_128BW_20000FS_200000NPKTS_100PLEN_10CR_3']
# test_list = ['SF_7N_128BW_20000FS_200000NPKTS_10PLEN_100CR_0']
test_list = ['SF_7N_128BW_20000FS_200000NPKTS_50PLEN_100CR_0']
# test_list = ['SF_7N_128BW_20000FS_200000NPKTS_10PLEN_100_doppler_apogee_simulation']
output_file_cnt = 0
SNR_dx = []
error_locs = []
for exp_folder in test_list:
    sf=7
    experiment_base_folder = exp_root_folder + '/' + exp_folder + '/' + 'error_out' + '*' + '.pkl'
    files_list = glob.glob(experiment_base_folder)
    print(len(files_list))
    
    output_file_cnt = 0
    for file in files_list:
        with open(file,'rb') as f: 
            #print(file)
            GND_TRUTH_PKT, OUTPUT_PKT, TOTAL_PKT_CNT, PKT_SNR = pickle.load(f)
            
            if len(OUTPUT_PKT) == len(GND_TRUTH_PKT): 
                #print(sum(np.subtract(GND_TRUTH_PKT,OUTPUT_PKT)))
                SF_x.append(sf)                
                print(int(np.count_nonzero(abs(np.subtract(GND_TRUTH_PKT,OUTPUT_PKT))) ))
                n_errors = int(np.count_nonzero(abs(np.subtract(GND_TRUTH_PKT,OUTPUT_PKT))) ) 
                if (n_errors > 0):
                    # plt.figure(1)
                    # plt.plot(GND_TRUTH_PKT)
                    # plt.plot(OUTPUT_PKT)
                    print("DX",np.argwhere(abs(np.subtract(GND_TRUTH_PKT,OUTPUT_PKT)) > 0))
                    error_locs.extend(np.argwhere(abs(np.subtract(GND_TRUTH_PKT,OUTPUT_PKT)) > 0))
                    plt.figure(2)
                    plt.plot(abs(np.subtract(GND_TRUTH_PKT,OUTPUT_PKT)))
                    plt.show()
                SER_Y.append(n_errors) 
                output_file_cnt += 1
                SNR_dx.append(PKT_SNR)
    SF_pkt_x.append(sf)
    PKT_CNT_Y.append(output_file_cnt) 

# plt.figure(5)
# plt.hist(error_locs)
# plt.show()
print('Packets Received:' ,PKT_CNT_Y)
# print(SNR_dx.shape)
nandx = np.argwhere(np.isnan(SNR_dx))
# if len(nandx) > 1:
    # print(nandx)
    # del SNR_dx[nandx]
    # del SER_Y[nandx]
    # del SF_x[nandx]


plt.figure(1)
plt.hist(SER_Y,bins=np.unique(SER_Y))
plt.xlabel('Number of Errors in a Packet')
plt.ylabel('Number Packets')

print("Mean SER:", np.mean(SER_Y))        
plt.figure(3)
plt.scatter(SNR_dx, SER_Y)
plt.xlabel('SNR')
plt.ylabel('Symbol Error Rate')
plt.title('Symbol Error Rate Vs. SNR')
plt.show()

# plt.figure(1)
# plt.scatter(SF_pkt_x, PKT_CNT_Y)
# plt.axhline(y=expected_pkt_cnt, color='r', linestyle='-')
# plt.title("Received Packet Count)")
# plt.xlabel("SF")
# plt.ylabel("# Received Packets")
# #plt.show()

# plt.figure(2)
# plt.scatter(np.arange(0, PKT_CNT_Y[0]), SER_Y)
# #plt.hist(SER_Y)
# plt.title("Symbol Errors")
# #plt.ylabel("freq")
# plt.ylabel("# Symbols != Ground Truth Packet")
# plt.xlabel("Received Pkt #")
# #plt.xlabel("# Symbols != Ground Truth Packet")
# plt.show()

print("Total symbol errrors: ", sum(SER_Y))