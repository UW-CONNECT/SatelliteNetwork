# Python_Realtime_TxRx
This is code for CSS-based modulation/demodulation using USRP B200 and GNURADIO. 

Realtime (~1s) Tx-Rx Mod/Demod Code with USRP + GNURADIO Python Implementation

## TODO 
* Add doppler correction code 
* Check different packet length sizes w.r.t. RAW_FS value in rx_stream
* Implement specific packet format including preamble, packet length, data, and end delimeter

## NOTES 
* Doppler correction can be done via GNURADIO using GPREDICT --> should we move this over? 
* Synthetic data generated are in simpleTX_sim/test_data/, and are named accordingly , make sure to change parameters in rx_stream.py

## Simple Test of Realtime Data Without USRPS (from synthetic data and GNURADIO channel model) 
1. simpleTX_sim/genTxFile.py generates raw CSS data into binaryTXdata according to the 'symbols' variable. 
2. Run simpleTX_fromUDP.grc 
3. Run rx_stream.py :: Packets should begin to be decoded in the output cmd line, and this should produce symbols + 1 

## Tx/Rx with Two USRPs instructions 
1. Connect each USRP to separate computers with GNURADIO and Python scripts. 
2. On one computer/USRP (1), run simpleRX_sim/simpleRx.grc in GNURADIO. Note the address and port 
   (127.0.0.1 and 4.9k): these should match the port on rx_stream.py. A QT GUI should pop up showing Rx noise. 
3. Run rx_stream.py on (1) and ensure that window shows "Packet not detected" if the transmitter is not on. 
4. On the other computer/USRP (2), run simpleTX_sim/simpleTX_fromUDP. Note the address and port (4444), 
   should also match the Python script. 
5. On (2), run the Python script simpleTX_sim/tx_stream_zmq.py: This transmits a 57-length packet 
   continuously, but this can be changed. 
6. Spectrogram on (1) should show a received signal and the window from rx_stream.py should show "Packet detected" 
   along with the demodulated symbols defined in tx_stream_zmq.py 

## simpleTx_sim Separate Instructions 
Run tx_stream_zmq.py and change port and ZMQ to match simpleTX_fromUDP.grc, as well as buffer size. 
Can view in GNURADIO GUI, save to file, or send over UDP again to rx_stream.py 

## simpleRx_sim Separate Instructions

### test_css_demod 
Standalone script for testing css_demod.py -> Just decodes first few symbols from the included data file. 

### rx_stream.py
This receives samples over UDP for CSS demodulation. This needs to be debugged and verified for 
realtime. 

Run simpleRX_sim_fromfile.grc - this repeats the waveform from rx_gen_sf7_small and sends it over UDP.
Then run rx_stream.py. 

Both test_css_demod.py and rx_stream.py use css_demod.py -> Just need to verify that packet detection
works across multiple queue 
