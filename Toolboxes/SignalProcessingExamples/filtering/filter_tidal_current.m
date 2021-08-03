close all; clear;

%% load tidal current data
%2008 11 24 00 05 00      0.76         261
%2008 11 24 00 15 00      0.32         262

load FPI0901.txt
data = FPI0901;

%% time
time = datenum(data(:,1),data(:,2),data(:,3),data(:,4),data(:,5),data(:,6));
dt_min = 10;
dt = 10*60;

%% current data
speed = data(:,7);
dir = data(:,8);
for i = 1:length(dir)
    if (dir(i) > 180)
        speed(i) = speed(i);
    else
        speed(i) = -speed(i);
    end
end

N = length(speed);

%% zero-padding
%speed=[speed' zeros(4*length(speed),1)'];
%N = length(speed);

%% plot
subplot(3,1,1)
grid
%plot(time-time(1),speed)
plot(speed)
%datetick('x',19)
%xlabel('Elapsed time (days)')

%% compute FFT
Y = fft(speed);
Y_mag=abs(Y);

subplot(3,1,2)
plot([1:N-1],Y_mag(2:end),'k')

%% setup frequency array
df = 1/(N*dt);
f0 = df;
fN = 1/(2*dt);
f = f0:df:fN;

subplot(3,1,3)
plot(f,Y_mag(2:N/2+1)/(N/2))
grid