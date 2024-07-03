#!/usr/bin/env python
'''
   Copyright 2015 Wolfgang Nagele

   Licensed under the Apache License, Version 2.0 (the "License");
   you may not use this file except in compliance with the License.
   You may obtain a copy of the License at

       http://www.apache.org/licenses/LICENSE-2.0

   Unless required by applicable law or agreed to in writing, software
   distributed under the License is distributed on an "AS IS" BASIS,
   WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
   See the License for the specific language governing permissions and
   limitations under the License.
'''
# from gnuradio import gr
import threading
import time
import socket


# class doppler_runner(threading.Thread):
  # def __init__(self, callback, gpredict_host, gpredict_port, verbose):
    # threading.Thread.__init__(self)

    # self.callback = callback
    # self.gpredict_host = gpredict_host
    # self.gpredict_port = gpredict_port
    # self.verbose = verbose


  # def run(self):
gpredict_host = 'localhost'
verbose = True 
gpredict_port = 7356
bind_to = (gpredict_host, gpredict_port)
server = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
server.bind(bind_to)
server.listen(0)

time.sleep(0.5) # TODO: Find better way to know if init is all done

while True:
  if verbose: print( "Waiting for connection on: %s:%d" % bind_to)
  sock, addr = server.accept()
  if verbose: print ("Connected from: %s:%d" % (addr[0], addr[1]))

  cur_freq = 0
  while True:
    data = sock.recv(1024)
    if not data:
      break

    if data[0] == 'F':
      print("Check new freq")
      freq = int(data[1:].strip())
      if cur_freq != freq:
        print("New frequency: %d" % freq)
        self.callback(freq)
        cur_freq = freq
      sock.sendall("RPRT 0\n")
    elif data[0] == 'f':
      sock.sendall("f: %d\n" % cur_freq)

  sock.close()
  if verbose: print( "Disconnected from: %s:%d" % (addr[0], addr[1]))

# if __name__ == "__main__":
    

# class doppler(gr.sync_block):
  # def __init__(self, callback, gpredict_host, gpredict_port, verbose):
    # gr.sync_block.__init__(self,
                           # name = "Gpredict Doppler",
                           # in_sig = None,
                           # out_sig = None)
    # doppler_runner(callback, gpredict_host, gpredict_port, verbose).start()