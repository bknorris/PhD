close all; clear;

% create ideal time array
dt_ideal = 1.0;
time_ideal = 0 : dt_ideal : 2000;

% create wave 1
T_wave1 = 20;
A_wave1 = 10;
wave1 = A_wave1*cos(2*pi*time_ideal/T_wave1);

% create wave 2
T_wave2 = 200;
A_wave2 = 1;
wave2 = A_wave2*cos(2*pi*time_ideal/T_wave2);

% combine the waves
wave_total_ideal = wave1 + wave2;

% plot it
plot(time_ideal,wave_total_ideal)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% OUR recording time,
mult = 21;
t = 0 : mult*dt_ideal : 2000;

% MY recorded wave
wave_recording = interp1(time_ideal,wave_total_ideal,t);

hold on
plot(t,wave_recording,'r.-')