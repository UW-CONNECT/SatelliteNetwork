'''
client :: Jackie Schellberg :: April 2024
Controls Tx of Chirps 
'''
# echo-client.py


import multiprocessing
from multiprocessing import Manager

from threading import Thread
from queue import Queue
import zmq
import socket
import threading 
import time
import numpy as np
import subprocess
# HOST = "127.0.0.1"  # Standard loopback interface address (localhost)

HOST = "10.139.37.113"  # The server's hostname or IP address
# HOST = "0.0.0.0"  # Standard loopback interface address (localhost)
PORT = 65432  # Port to listen on (non-privileged ports are > 1023)

HEADER = 64
# SERVER = ""
# Another way to get the local IP address automatically
SERVER = socket.gethostbyname(socket.gethostname())
ADDR = (SERVER, PORT)
FORMAT = 'UTF-8'
DISCONNECT_MESSAGE = "!DISCONNECT"
conns = []
addrs = []
DOPPLER_PROC_TIME = 20 * 60 * 2 # time to go through the stored doppler file. double check this 
EXP_ROOT_FOLDER = '../../experiment_data_synced'
#exp_root_folder = '../../simulated_data'
#exp_folder = str(NUMPKTS) + '_pkts_' + str(SF) + '_SF'  
# exp_folder = 'SF_' + str(SF) + 'N_' + str(N) + 'BW_' + str(BW) + 'FS_' + str(FS) +\
# 'NPKTS_' +str(NUMPKTS) + 'PLEN_' +st+ str(N) + 'BW_' + str(BW) + 'FS_' + str(FS) +\
# 'NPKTS_' +str(NUMPKTS) + 'PLEN_' +str(PAYLOAD_LEN) +'CR_' +str(CR)
# exp_folder = str(trial_num_name_out)
RAW_DATA_FILE_NAME = "trial1"

START_SYMBOL = str(1) # to trigger the start of data collection 
END_SYMBOL = str(0)   # to trigger early exit from data collection 
def server_action(conn, addr):
    print(f"[NEW CONNECTION] {addr} connected.")
    connected = True
    while connected:
        msg_length = conn.recv(HEADER).decode(FORMAT)
        if msg_length:
            msg_length = int(msg_length)
            msg = conn.recv(msg_length).decode(FORMAT)
            if msg == DISCONNECT_MESSAGE:
                connected = False
            print(f"[{addr}] {msg}")
        conn.send("Msg received".encode(FORMAT))

    conn.close()

def send_tx_zmq_data(trial_number): 
    '''
    Send doppler data over ZMQ to the doppler transponder. 
    '''
    FILE_PATH = EXP_ROOT_FOLDER + '\\' + str(trial_number) + '\\' + RAW_DATA_FILE_NAME
    context = zmq.Context() 
    # zsocket = context.socket(zmq.PUB) 
    zsocket = context.socket(zmq.PUB) 
    zsocket.bind("tcp://127.0.0.1:4444")
    experiment_start_time = time.time() # [TODO] save this to synchronize the experiments
    print("Started sending TX data over ZMQ. Start time: ", experiment_start_time)
    with open(FILE_PATH, "rb") as f:
        # nbytes= int((32 /8) *2 * 200000)
        # print(nbytes)
        nbytes = 200000
        while (byte := f.read(nbytes)):
            # Do stuff with byte.        
            
            # byte = np.frombuffer(byte, dtype=np.float32)
            zsocket.send(byte)
            print("Sent byte")
            # time.sleep(.1)
    print("Done sending the data over ZMQ")

def run_sync_experiment(): 
    '''
    Runs experiments within the time of the doppler curve
    '''
    current_trial_number = 0
    start_time=time.time()
    while(time.time() - start_time < DOPPLER_PROC_TIME): 
    # send multiple packets across the doppler curve, make sure they wait long enough 
        proc_file_name = EXP_ROOT_FOLDER + '/' + str(current_trial_number) + '/' + 'proc_file_time'
        f = open(proc_file_name)
        DATA_FILE_PROC_TIME = f.read()
        f.close()
        print("Current time to run for this one: ", DATA_FILE_PROC_TIME)
        time.sleep(8)  # give time for python on Rx to startup
        
        send_tx_zmq_data(current_trial_number)
        current_trial_number = current_trial_number + 1 
        time.sleep(float(DATA_FILE_PROC_TIME)+2)

if __name__ == "__main__":    
    status = [END_SYMBOL ]
    
    # procHandle = run(DOPPLER_GNURADIO)
    # os.system("python3 " + DOPPLER_GNURADIO +"&") 
    # time.sleep(10) # give GNURADIO some time to wakeup
    with socket.socket(socket.AF_INET, socket.SOCK_STREAM) as s:
        print("A socket.")
        s.connect((HOST, PORT))
        s.sendall(b"RECEIVER RX CLIENT CONNECTED!!!")
        
        recver = list()
        counter = 0

        while True:
            # if s.poll(10) != 0:
            message = s.recv(1024).decode()  
            print('Received from server: ' + message[0])
            if (message[0] == START_SYMBOL): 
                print("Bouta start")
                status = [START_SYMBOL ]
                thread = threading.Thread(target=run_sync_experiment, args=())
                thread.start()
            elif (message[0] == END_SYMBOL): 
                status = [END_SYMBOL ]
        
    print(f"Received {data!r}")
