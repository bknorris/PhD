% Script used to do the real-time spectrogram demo during the Short Course 
% Friday session. Running this script grabs the handle to the default
% microphone and records 4[s] of audio. 

close all
clear all
%
N = 2048*1; Ninc = round(N/8);  %  overlap (Doppler)
% N = 2048*8; Ninc = round(N/8);  %  overlap (Fan)
%
nbits = 16; nchan = 1; fs = 12000;
Trec = 5;  inputID = -1;
%
rr = audiorecorder(fs, nbits, nchan, inputID);
%
input('Press enter to start recording')
recordblocking(rr, Trec)
pause(0.1)
yy = getaudiodata(rr);
%
dt = 1/fs; T = N*dt; df = 1/T;
%
%  Hanning window
[W, Aw] = hanning_window(N);
W = W.';
%
yy2D = fill_2D_array(yy, N, W, 8);
psd2D = Aw*Gxx_2D(yy2D, fs);
%
[Npsd, Nrec] = size(psd2D);
Nrange = 2:round(0.8*Npsd);   %  Select frequency range to display
psd2D = 10*log10(psd2D(Nrange,:));
%
dt_overlap = Ninc*dt;
times = ((1:Nrec) - 0.5)*dt_overlap; % Make times the middle of the records
freqs = (Nrange - 1)*df;
%
position = [50, 50, 800, 600];
figure(1)
set(gcf, 'position', position)
hh = plot(yy);
%
figure(2)
set(gcf, 'position', position)
%
% upper_dB = -20;  lower_dB = -100;
% psd2D(psd2D > upper_dB) = upper_dB;
% psd2D(psd2D < lower_dB) = lower_dB;
%
% % hh = spectroplot(times, freqs, psd2D);
hh = imagesc(times, freqs, psd2D);
axis xy
colormap(jet)
colorbar
set(gca, 'fontsize', 18)
xlabel('Time [s]', 'fontsize', 16)
ylabel('Frequency [Hz]', 'fontsize', 16)
%
figure(3)
set(gcf, 'position', position)
%
% upper_dB = -20;  lower_dB = -100;
% psd2D(psd2D > upper_dB) = upper_dB;
% psd2D(psd2D < lower_dB) = lower_dB;
%
% % hh = spectroplot(times, freqs, psd2D);
hh = imagesc(times, freqs, psd2D);
axis xy
colormap(jet)
colorbar
set(gca, 'fontsize', 18)
xlabel('Time [s]', 'fontsize', 16)
ylabel('Frequency [Hz]', 'fontsize', 16)
axis([0 Trec 800 1200])
%
