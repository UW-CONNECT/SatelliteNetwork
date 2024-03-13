% file_dir = 'C:\\Users\\schellberg\\Downloads\\USRP_transponder_SNR_test';
close all; clear all; clc;
addpath('C:\Users\schellberg\Documents\schellberg\tutorials\USRP\LoRa_Learn_Simple')
% file_dir = 'single10_SF9'
%file_dir = '../../simulated_data_SNR/sf7_noDop_SNR1-5hundredth_4dB';
file_dir = '../../experiment_data/snr_xlater_filter/snr_2ndfloor_SF7';
x = load_file(file_dir,1,0);
st_dx = 50000;
end_dx = st_dx+ 300000;
% end_dx= 13635600;
x = x(st_dx:end_dx);
figure(1);
% plot(10*log10(abs(x)))
fs = 200000;
x = lowpass(x, 20e3, fs);
plot(10*log10(abs(x)))
% plot(mag2db(abs(x)))
% pspectrum(x(185483:end),'spectrogram')
%plot(abs(x))
% pspectrum(x(1:300000),'spectrogram')
% x = lowpass(x, 25e3);
% start_pow = 180000
start_pow = 290000/2 + 1
hold on;
xline(start_pow)

nsamps = length(x(start_pow:end));

xline(nsamps,'r')
% np = sum(10*log10(abs(x(1:nsamps))));
% sp = sum(10*log10(abs(x(185483:end))));
% np = sum(abs(x(1:nsamps)))/nsamps;
% sp = sum(abs(x(185483:end)))/nsamps;
% x = x/max(x)
np = rms(x(1:nsamps))^2; %/nsamps 
% np2 = bandpower(x,fs,[0 fs/2])
sp = rms(x(start_pow:end))^2; %/nsamps 
% sp = sum(x())
% SNR = 10*log(sp - np / np) 
SNR = 10*log((sp-np) / np) 