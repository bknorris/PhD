function [time_clip,fs_clip_detrend]=get_wake_data(raw_data);

%load 030514_Testing.txt

%% calibration:  BLUE gauge
depth = [10 20 30 40];
ideal_count_array=0:1:4000;
blue_cal_counts=[804 1576 2453 3083];
P_blue=polyfit(blue_cal_counts,depth,1);
cal_curve_blue = polyval(P_blue,ideal_count_array);

%% convert counts to free surface
%raw_data = X030514_Testing';

[m n] = size(raw_data);

counts = reshape(raw_data',m*n,1);
free_surface=P_blue(1)*counts + P_blue(2);

%% make time array:  Y14,M03,D05,H07,M11,S00
start_time = datenum(2014,3,5,7,11,00);
sampling_f = 10;
dt_sec = 1/sampling_f;
dt = dt_sec/60/60/24;
time = start_time:dt:start_time+dt*(m*n-1);

%% clip and detrend
fs_clip_detrend=detrend(free_surface(40300:50850));
fs_clip_detrend(4001:end)=[];
time_clip = time(40300:50850);
time_clip(4001:end)=[];
N = length(time_clip);
plot(fs_clip_detrend)
grid

