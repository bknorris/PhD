close all; clear;

%dt = 0.37;
%t = dt:dt:12;
t = linspace(0,36,2^8);
dt = t(2)-t(1);
%dt = 0.36;
%t = 0:dt:12;
T1 = 12;
T2 = 4;
T3 = 1.6;
omega1 = 2*pi/T1;
omega2 = 2*pi/T2;
omega3 = 2*pi/T3;
amp1 = 1.6;
amp2 = 0.85;
amp3 = 0.4;
eta = amp1*sin(omega1*t) + amp2*cos(omega2*t) + amp3*cos(omega3*t);% + rand(1,length(t));
%subplot(2,1,1)
plot(t,eta)
N = length(eta);

%% get Fourier coefficents
pause
tstart=tic;
[A0,A,B] = get_fourier_coeffs(eta,1);

%% get energy from Fourier coefficents
[S,amp_sq,f]=get_energy(A,B,N,dt);
time_fourier = toc(tstart);

%% FFT
tstart=tic;
Xn = fft(dt*eta,N);
%amp_sq_fft = 1/(dt*N) * (Xn.*conj(Xn));
%amp_sq_fft(N/2+2:end) = [];
Gyy = 2/(dt*N) * (Xn.*conj(Xn));
Gyy(N/2+2:end) = [];
time_fft=toc(tstart);

%% plot and compare traditional Fourier to FFT Amplitudes
% figure
% subplot(2,1,1)
% bar(f,amp_sq,'r')
% grid
% title('Traditional Fourier')
% ylabel('Amplitude Squared (units^2)')
% 
% subplot(2,1,2)
% bar(f,amp_sq_fft(2:N/2+1),'k')
% grid
% xlabel('f (Hz)')
% title('FFT')
% ylabel('Amplitude Squared (units^2)')

%% plot and compare traditional Fourier to FFT Power spectra
pause
figure
%subplot(2,1,1)
plot(f,S,'r');
grid
ylabel('Energy Density (units^2/Hz or units^2*s)')
title('One-sided power spectrum:  Traditional Fourier')
%disp(['Total signal energy is ',num2str(trapz(f,S)), ' units^2'])

% subplot(2,1,2)
% plot(f,Gyy(2:N/2+1),'k')
% grid on
% xlabel('f (Hz)')
% ylabel('Energy Density (units^2/Hz or units^2*s)')
% title('One-sided power spectrum:  FFT version')

%% speed-up
disp(['FFT is ',num2str(time_fourier/time_fft),' times faster than traditional Fourier.'])
