close all; clear;

load X030514_Testing.txt;

%% get wake data
[time,eta]=get_wake_data(X030514_Testing);
N = length(eta);
% plot(time,eta)
% grid
% datetick('x',13)
% ylabel('Free surface (cm)')
% xlabel('Time')

%% hanning window
h = hanning(N);
eta = eta .* h;

%% zero pad
eta = [eta; zeros(100000,1)];
N = length(eta);

%% FFT 
Y = fft(eta);
Y_mag=abs(Y)/(N/2);
Y_one_side=Y_mag(2:N/2+1);

fs = 10;
dt = 0.1;
f0 = fs/N;
df = f0;
fN = fs/2;
f = f0:df:fN;

figure
plot(f,Y_one_side)
grid