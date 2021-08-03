% This tutorial attempts to explain how zero-padding works.
% You should already have some exposure to how the DFT works and
% be have seen the effects of zero paddding and windowing.
%
close all; clear;

fs = 100;
dt = 1/fs;
N = 2^8; % Create a signal with N samples = 1 second
N_orig = N;
n = 0:N-1; % sample numbers
T = (N-1)*dt; % sampling times over a one second duration.
t = 0:dt:T;
f1 = 3.0;
x = 1.0*cos(2*pi*f1*t/T);

%% windowed x
%h=hanning(N);
%x=x.*h';

%% padding x with 1000 zeros
%x = [x, zeros(1,1000)];
%N = length(x);
%n =0:N-1; 

%plot the magnitude spectrum 
X_mags = abs(fft(x));

%% generate frequency array
%f0 = 1/(N*dt);
%fN = 1/(2*dt);
%df = f0;
%f = f0:df:fN;
%f = [0 f]; 

subplot(2,1,1)
plot(n,x)
xlabel(['Samples,n of N = ',num2str(N)]);
ylabel('Amplitude')
title('Time-Domain Signal');
 
subplot(2,1,2)
plot(n(1:N/2+1),X_mags(1:N/2+1)/(N_orig/2))
xlabel('Frequency bin index, k');
ylabel('Amplitude (units)');
title('Frequency Content Magnitudes:  one-sided');

% return
% 
% num_disp_bins = 30;
% subplot(2,1,2)
% plot([0:num_disp_bins-1], X_mags(1:num_disp_bins));
% hold on
% plot([0:num_disp_bins-1], X_mags(1:num_disp_bins),'k.');
% hold off
% xlabel('Frequency Bins');
% ylabel('Magnitude');
% title('Frequency Content Magnitudes');

%Illustrate Basis functions
figure
plot(n,x,'k')
xlabel('Samples');
ylabel('Amplitude')
hold on
for k = 0 : 70 
    plot(x,'kx')
    hold on
    cos_basis = cos(2*pi*k*n/N);
    sin_basis = sin(2*pi*k*n/N);
     
    plot(n,cos_basis,'r')
    plot(n,sin_basis,'g')
    title({['Signal being anlaysed (black) with basis functions'...
        ' that have ' num2str(k) ' cyles over '  ...
        num2str(N) ' samples']...
        ['Cosine basis function shown in red, Sine in green']})
    hold off
    pause
end