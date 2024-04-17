'''
For synchronzing experiments: 
https://github.com/MattCrook/python_sockets_multi_threading/blob/master/server.py
https://realpython.com/python-sockets/
https://docs.python.org/3.5/library/subprocess.html#subprocess.Popen
'''
# echo-server.py

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
HOST = "0.0.0.0"  # Standard loopback interface address (localhost)
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
 
def client_action(conn, addr):
    # with conn:
        # print(f"Connected by {addr}")
        # while True:
            # data = conn.recv(1024)
            # print("Waiting")
            # if not data:
                # break
            # conn.sendall(data)

    print(f"[NEW CONNECTION] {addr} connected.")
    connected = True
    while connected:
        # msg_length = conn.recv(HEADER).decode(FORMAT)
        msg_length = conn.recv(1024).decode()
        msg = msg_length
        msg_length = len(msg)
        if msg_length:
            # msg_length = int(msg_length)
            # msg = conn.recv(msg_length).decode(FORMAT)
            # if msg == DISCONNECT_MESSAGE:
                # connected = False
            print(f"[{addr}] {msg}")
        conn.send("Msg received".encode(FORMAT))

    conn.close()

def send_tx_zmq_data(trial_number): 
    '''
    Send doppler data over ZMQ to the doppler transponder. 
    '''
    FILE_PATH = EXP_ROOT_FOLDER + '\\' + str(trial_number) + '\\' + RAW_DATA_FILE_NAME
    context = zmq.Context() 
    zsocket = context.socket(zmq.PUB) 
    zsocket.bind("tcp://127.0.0.1:4444")
    experiment_start_time = time.time()
    print("Started sending TX data over TCP. Start time: ", experiment_start_time)
    with open(FILE_PATH, "rb") as f:
        nbytes= int((32 /8) *2 * 200000)
        print(nbytes)
        while (byte := f.read(nbytes)):
            # Do stuff with byte.        
            zsocket.send(byte)
            byte = np.frombuffer(byte, dtype=np.float32)
            #print("Sent byte")
            time.sleep(.1)
    print("Done sending the data over ZMQ")

def run_experiment(conns,addrs): 
    '''
    Tells when to start experiment to the clients. 
    '''
    
    # check how many experiment files we need to check -- all in one shared directory, named based on trial #   
    current_trial_number = 0
    while(True):
        #print(len(conns))
        if (len(conns) >= 2):
            start_time = time.time()    # Using POSIX time to syncronize the doppler curves after the fact. 
            print("START RUNNING EXPERIMENT:")        
        # 
            experiment_start = '1_experiment'
            time.sleep(1)
            for c in range(0, len(conns)):
                conn = conns[c]
                conn.send(experiment_start.encode(FORMAT))
                print("Tried to send a message.")
            time.sleep(1)
            #break;
            # os.system("python3 " + DOPPLER_GNURADIO +"&") 
            while(time.time() - start_time < DOPPLER_PROC_TIME): 
                proc_file_name = EXP_ROOT_FOLDER + '/' + str(current_trial_number) + '/' + 'proc_file_time'
                f = open(proc_file_name)
                DATA_FILE_PROC_TIME = f.read()
                f.close()
                
                # Run our streamer with the correct configuration. 
                p = subprocess.Popen("python rx_stream_from_args.py " + str(current_trial_number) )
                time.sleep(float(DATA_FILE_PROC_TIME)+2)
                p.kill()
                
                print("Killed this process.")
                current_trial_number = current_trial_number + 1 
            # while(time.time() - start_time < DOPPLER_PROC_TIME): 
                # # send multiple packets across the doppler curve, make sure they wait long enough 
                # proc_file_name = EXP_ROOT_FOLDER + '/' + str(current_trial_number) + '/' + 'proc_file_time'
                # f = open(proc_file_name)
                # DATA_FILE_PROC_TIME = f.read()
                # f.close()
                # print("Current time to run for this one: ", DATA_FILE_PROC_TIME)
                # send_tx_zmq_data(current_trial_number)
                # current_trial_number = current_trial_number + 1 
                # time.sleep(DATA_FILE_PROC_TIME+2)
if __name__ == "__main__":
    thread = threading.Thread(target=run_experiment, args=(conns, addrs))
    thread.start()
    with socket.socket(socket.AF_INET, socket.SOCK_STREAM) as s:
        s.bind((HOST, PORT))
        s.listen()
        print(f"[LISTENING] Server is listening on {SERVER}")
        while True:
            conn, addr = s.accept()
            conns.append(conn)
            addrs.append(addr)
            thread = threading.Thread(target=client_action, args=(conn, addr))
            thread.start()
            print(f"[ACTIVE CONNECTIONS] {threading.activeCount() - 1}")