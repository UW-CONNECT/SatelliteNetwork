'''
transponder_client :: Jackie Schellberg :: April 2024 

Transponder client. Just listens for a message from the server to trigger start/stop 
of the doppler transponder. GNURADIO is expected to already be running; ZMQ connection is triggered from here. 

'''
# echo-client.py

import socket
import numpy as np
import subprocess
import shlex
import zmq
import os 
import time
import threading 
# HOST = "127.0.0.1"  # The server's hostname or IP address
HOST = "10.139.37.113"  # The server's hostname or IP address
# HOST = "192.168.56.1"
PORT = 65432  # The port used by the server
FORMAT = 'UTF-8'
DOPPLER_GNURADIO = "doppler_transponder_noGUI.py"  # name of script to run (generated from GNURADIO)
# FILE_PATH = r'J:\schellberg\indoor_exp_feb_2024\doppler_simulation_files\GNURADIO_linear_8k_sep'
# FILE_PATH = r'..\..\indoor_exp_feb_2024\doppler_simulation_files\GNURADIO_linear_8k_sep'
# FILE_PATH = r'../../indoor_exp_feb_2024/doppler_simulation_files/GNURADIO_linear_8k_sep'
FILE_PATH = r'..\..\doppler_simulation_files\GNURADIO_linear_8k_sep'
START_SYMBOL = str(1) # to trigger the start of data collection 
END_SYMBOL = str(0)   # to trigger early exit from data collection 
# status = [END_SYMBOL ]
def server_action(conn, addr):
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
        msg_length = conn.recv(HEADER).decode(FORMAT)
        if msg_length:
            msg_length = int(msg_length)
            msg = conn.recv(msg_length).decode(FORMAT)
            if msg == DISCONNECT_MESSAGE:
                connected = False
            print(f"[{addr}] {msg}")
        conn.send("Msg received".encode(FORMAT))

    conn.close()
    
def run(script):
    scriptArgs = shlex.split(script)
    commandArgs = ["python3"]
    commandArgs.extend(scriptArgs)
    procHandle = subprocess.Popen(commandArgs, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
    return procHandle
    
def start_GNURADIO(status):
    print("Controlling doppler.")
    
    
    # procHandle = None
    # while (True):
        # print(status[0])
        # if (status[0] == START_SYMBOL):# and procHandle.poll() is None):
            # procHandle = run(DOPPLER_GNURADIO)
            # print("Started======")
            # time.sleep(5)
            # isScriptRunning(procHandle)
            # stopScript(procHandle)
            # print getOutput(procHandle)
            # status = [str(2)] # set process only once 
        # elif (procHandle is not None and status == END_SYMBOL): 
            # procHandle.terminate()
        # time.sleep(.1)

def send_doppler_socket():
    '''
    Send doppler data over ZMQ to the doppler transponder. 
    '''
    context = zmq.Context() 
    zsocket = context.socket(zmq.PUB) 
    zsocket.bind("tcp://127.0.0.1:4444")
    
    # save the start time for alignment 
    experiment_start_time = time.time()
    doppler_start_file = "doppler_start_file"
    f = open(doppler_start_file, "w")
    f.write(str(experiment_start_time))
    f.close()
    
    print("Started sending doppler data over TCP. Start time: ", experiment_start_time)
    with open(FILE_PATH, "rb") as f:
        # nbytes= int((32 /8) *2 * 200000)
        nbytes= 200000
        print(nbytes)
        while (byte := f.read(nbytes)):
            # Do stuff with byte.        
            zsocket.send(byte)
            byte = np.frombuffer(byte, dtype=np.float32)
            #print("Sent byte")
            time.sleep(.1)
    print("Done sending the data over ZMQ")
    # pass      
if __name__ == "__main__":    
    status = [END_SYMBOL ]
    
    # procHandle = run(DOPPLER_GNURADIO)
    # os.system("python3 " + DOPPLER_GNURADIO +"&") 
    # time.sleep(10) # give GNURADIO some time to wakeup
    with socket.socket(socket.AF_INET, socket.SOCK_STREAM) as s:
        print("A socket.")
        s.connect((HOST, PORT))
        s.sendall(b"TRANSPONDER CONNECTED!!!! ")
        
        recver = list()
        counter = 0

        while True:
            # if s.poll(10) != 0:
            message = s.recv(1024).decode()  
            print('Received from server: ' + message[0])
            if (message[0] == START_SYMBOL): 
                print("Bouta start")
                status = [START_SYMBOL ]
                thread = threading.Thread(target=send_doppler_socket, args=())
                thread.start()
            elif (message[0] == END_SYMBOL): 
                status = [END_SYMBOL ]
        
    print(f"Received {data!r}")