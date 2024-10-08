'''
For processing the experiment data into plots and such, given ground truth data. 
'''
import pickle 
import glob 
import numpy as np
from matplotlib import pyplot as plt
SF = [7, 8, 9, 10, 11, 12]
BW = 20000
FS = 200000
NUMPKTS = 10 
PAYLOAD_LEN = 100
#exp_root_folder = 'indoor_testing'
# exp_root_folder = '../experiment_data_maxtx1transponder'
#exp_root_folder = '../simulated_data'
#exp_root_folder = '../simulated_data'
exp_root_folder = '../experiment_data'
SF_x = []
SNR = []
SF_pkt_x = []
PKT_CNT_Y = []
SER_Y = []
expected_pkt_cnt = 10 # this should be saved in our ground truth file... 
#for sf in SF: 
for sf in [7]:
    #exp_folder = '10_pkts_'+str(sf)+'_SF'
    N = 2**sf 
    exp_folder =  'SF_' + str(sf) + 'N_' + str(N) + 'BW_' + str(BW) + 'FS_' + str(FS) +\
'NPKTS_' +str(NUMPKTS) + 'PLEN_' +str(PAYLOAD_LEN) + 'CR_0'+'*'
    trial_name = "trial1"
    experiment_base_folder  = exp_root_folder + '/' + exp_folder + '/' + 'error_out' + '*' + '.pkl'
    print(experiment_base_folder)
    files_list = glob.glob(experiment_base_folder)
    print(len(files_list))
    
    output_file_cnt = 0
    for file in files_list:
        with open(file,'rb') as f: 
            print(file)
            GND_TRUTH_PKT,  OUTPUT_PKT, TOTAL_PKT_CNT, PKT_SNR = pickle.load(f)
            
            if len(OUTPUT_PKT) == len(GND_TRUTH_PKT): 
                #print(sum(np.subtract(GND_TRUTH_PKT,OUTPUT_PKT)))
                SF_x.append(sf)   
                SNR.append(PKT_SNR)
                print(int(np.count_nonzero(abs(np.subtract(GND_TRUTH_PKT,OUTPUT_PKT))) ))
                SER_Y.append(int(np.count_nonzero(abs(np.subtract(GND_TRUTH_PKT,OUTPUT_PKT))) )) 
                output_file_cnt += 1
    SF_pkt_x.append(sf)
    PKT_CNT_Y.append(output_file_cnt) 

print("Mean Symbol Errors:",np.mean(SER_Y))

plt.figure(1)
plt.hist(SER_Y)
plt.xlabel('Number of incorrect symbols')
plt.ylabel('Frequency')
ax = plt.gca()
# ax.set_xlim([-20, 20])
ax.set_ylim([0, 30])

plt.figure(2)
plt.scatter(SNR, SER_Y)
ax = plt.gca()
ax.set_xlim([-20, 20])
ax.set_ylim([0, 110])
plt.xlabel('SNR (dB)')
plt.ylabel('# Symbol Errors')
plt.show()

# plt.figure(1)
# plt.scatter(SF_pkt_x, PKT_CNT_Y)
# plt.axhline(y=expected_pkt_cnt, color='r', linestyle='-')
# plt.title("Received Packet Count)")
# plt.xlabel("SF")
# plt.ylabel("# Received Packets")
# #plt.show()

# plt.figure(2)
# plt.scatter(SF_x, SER_Y)
# #plt.hist(SER_Y)
# plt.title("Symbol Errors")
# #plt.ylabel("freq")
# plt.ylabel("# Symbols != Ground Truth Packet")
# plt.xlabel("SF")
# #plt.xlabel("# Symbols != Ground Truth Packet")
# plt.show()