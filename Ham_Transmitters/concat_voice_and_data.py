import numpy as np 


outfile_name = 'KJ5DKX_SF7ChirpTrain_5smallpkts'
voice_file = 'KJ5DKX_SSB_Modulated_voice_200k_Experimental'   
# chirp_file = '../../HAM_TESTS/0HzS_SF_7N_128BW_2500FS_200000NPKTS_1PLEN_100CR_0/trial1' 
chirp_file = '../../HAM_ONES_TEST/0HzS_SF_7N_128BW_2500FS_200000NPKTS_5PLEN_10CR_0/trial1' 

v1 = open(voice_file, 'rb')
v1_x = v1.read()

c1 = open(chirp_file, 'rb')
c1_x = c1.read()

print(len(v1_x))
print(len(c1_x))

# outdata = np.concatenate(( v1_x, c1_x ))
outdata= v1_x + c1_x

file = open(outfile_name, 'wb')
file.write(outdata)
file.close()

# file = open(outfile_name, 'wb')
d1 = open(outfile_name, 'rb')
d1_x = d1.read()

print(len(d1_x))