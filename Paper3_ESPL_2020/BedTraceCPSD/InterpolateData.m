%Interpolate VP data and Aquadopp data for spectral analysis.
%
% Paper 3: Sediment Motion in Mangroves
%
% This is Version 1.0 of this script
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear
close all
%Load Files
expname = 'F2F3_1';
dir1 = 'd:\Projects\Mekong_W2015\DataAnalysis\Paper3\VPs\';
files = {'11March2015_Vels.mat';'11March2015_Sen.mat'};
fn1 = whos('-file',[dir1 files{1}]);fn1 = {fn1.name};
V = matfile([dir1 files{1}]);
dir2 = 'd:\Projects\Mekong_W2015\DataAnalysis\Paper3\BottomTrackCPSD\';
dir3 = 'd:\Projects\Mekong_W2015\Data\Aquadopp\F2F3\';
bdfiles = dir([dir2 '*bdtrace.mat']);bdfiles = {bdfiles.name};
aqfiles = dir([dir3 '*2015.mat']);aqfiles = {aqfiles.name};
aqorder = [3 1 2];
%%%
start = datenum('11-Mar-2015 15:05:05');
stop = datenum('11-Mar-2015 16:50:11');
%%%

for ii = 1:3
    %Load VP data
    load([dir2 bdfiles{ii}]);
    vp = V.(fn1{ii});
    vpt = vp.time;id = find(vpt>=start*vpt<=stop);vpt = vpt(id);
    x = nanmean(vp.x(id,13:17),2); %cross-shore
    y = nanmean(vp.y(id,13:17),2); %along-shore
    
    %Load Aquadopp data, compute depth from P
    load([dir3 aqfiles{aqorder(ii)}])
    ll = (sin(aqdp.metadata.lat/57.29578))^2;
    at = aqdp.datenum;id = find(at>=start&at<=stop);at = at(id);
    P = aqdp.pressure(id);
    P = cmgbridge(P,100,100,1000);
    P(isnan(P)) = 0;
    g = zeros(length(P),1);h = zeros(length(P),1);
    for j = 1:length(P)
        g(j,:) = 9.780318*(1.0+(5.2788E-3+2.36E-5*ll)*ll)+1.092E-6*P(j,:);
        h(j,:) = ((((-1.82E-15*P(j,:)+2.279E-10)*P(j,:)-2.2512E-5)*P(j,:)+9.72659)*P(j,:))/g(j,:);
    end
    
    %create new timebase for all measurements, up/downsample to 25 Hz
    tb = linspace(start,stop,floor(length(x)/2)); %25Hz time series
    xi = linspace(1,length(h),floor(length(x)/2));
    depth = interp1(1:length(h),h,xi,'linear');
    xi = linspace(1,length(bd),floor(length(x)/2));
    bd25 = interp1(1:length(bd),bd,xi,'linear');
    xi = linspace(1,length(x),floor(length(x)/2));
    x = interp1(1:length(x),x,xi,'linear');
    xi = linspace(1,length(y),floor(length(y)/2));
    y = interp1(1:length(y),y,xi,'linear'); 
    
    %preprocess, remove nans
    tb(isnan(bd25)) = [];
    depth(isnan(bd25)) = [];
    x(isnan(bd25)) = [];
    y(isnan(bd25)) = [];
    bd25(isnan(bd25)) = [];
    
    %force inputs to be even
    if mod(length(x),2) == 1
        x = x(1:length(x)-1);
    end
    if mod(length(y),2) == 1
        y = y(1:length(y)-1);
    end
    if mod(length(depth),2) == 1
        depth = depth(1:length(depth)-1);
    end
    if mod(length(bd25),2) == 1
        bd25 = bd25(1:length(bd25)-1);
    end
    
    %Set up spectra settings
    fs = 25;
    win = 1200; %seconds
    step = 60; %seconds
    nwin = fs*win;
    avt = fs*step;
    nsamp = length(tb);
    
    %load data to run processing on into structure
    data.h = depth;
    data.bd = bd25;
    data.x = x;
    data.y = y;
    w = 2000;

    %run spectra
    vars = compute_spectra(tb,data,win,step,w,fs);
    
    %plot
    fname = [expname '_' fn1{ii}];
    plot_spec_results(fname,vars)
end
