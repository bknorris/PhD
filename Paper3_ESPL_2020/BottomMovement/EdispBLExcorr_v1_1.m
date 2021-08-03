%Plot TKE & BLE (bed movement) for an experiment. Then, calculate
%xcorr(TKE,BLE) if sensible.
%
% Paper 3: Sediment Motion in Mangroves
%
% This is Version 1.1 of this script
%
% V 1.0: After attempting multiple methods to generate cross-correlations,
% I have landed on a method. This code produces estimates of the time-lag
% between TKE dissipation and BLE, and additionally computes a scatter plot
% of bed movement, dissipation, and bed variance.
% V 1.1: Unable to produce reliable correlations, I am going to employ a
% method for resolving time series using MA, AM or ARIMA processes on
% non-stationary t-s (bdh). Will compare dissipation, wave energy vs. bdh.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear
close all
expname = 'VTA';
start = datenum('13-Mar-2015 00:00:00');
stop = datenum('15-Mar-2015 00:00:00');
wbbl = [0.0025 0.0025]; %from WaveBoundaryThickness.m
bins = 13:17;
hoffset = 99; %difference between vp heading and transect (deg)
%
dir1 = 'd:\Projects\Mekong_W2015\DataAnalysis\Paper3\BottomTrack\';
bdfiles = {'VTA1a_bdtrace.mat'};
dir2 = 'd:\Projects\Mekong_W2015\DataAnalysis\Paper3\VPs\';
vfiles = {'13March2015a_Vels.mat'};
dir3 = 'd:\Projects\Mekong_W2015\DataAnalysis\Paper3\Turbulence\';
tfiles = {'13March2015a_VelsTKE.mat'};
%
%Load Wave Energy files
load('d:\Projects\Mekong_W2015\DataAnalysis\Paper3\WaveEnergy\VTA1a_duetwvs.mat');
%
%Plot colors
if strfind(expname,'VTA')
    cl = [0 0 0];
else
    cl = brewermap(5,'Blues');
    cl = cl(3:5,:);
end
%Rangebins height
rb = [0.0400,0.0410,0.0419,0.0429,0.0439,...
    0.0448,0.0458,0.0467,0.0477,0.0487,...
    0.0496,0.0506,0.0516,0.0525,0.0535,...
    0.0545,0.0554,0.0564,0.0574,0.0583,...
    0.0593,0.0602,0.0612,0.0622,0.0631,...
    0.0641,0.0651,0.0660,0.0670,0.0680,...
    0.0689,0.0699,0.0708,0.0718,0.0728];
sdir = 'd:\Projects\Mekong_W2015\Figures\Paper3\BottomTracking\';

%Hone in on a specific time period when TKE is high; look for xcorrelation
%between the two signals
t1 = {'13-Mar-2015 15:00:00'};
t2 = {'13-Mar-2015 19:00:00'};
for i = 1:length(bdfiles)
    start = datenum(t1{i});
    stop = datenum(t2{i});
    %Load the data
    bd = load([dir1 bdfiles{i}]);
    fn = whos('-file',[dir2 vfiles{i}]);fn = {fn.name};
    v = matfile([dir2 vfiles{i}]);
    load([dir3 tfiles{i}]); %Turbulence
    if strfind(expname,'VTA')
        ni = 1;
    else
        ni = 1:3;
    end
    for ii = ni
        %Prepare bottom tracking data
        cp1 = find(bd.(fn{ii}).time >= start & bd.(fn{ii}).time <= stop);
        bdt = bd.(fn{ii}).time(cp1);
        bdh = bd.(fn{ii}).bdist(cp1);
        bid = find(~isnan(bdh),1,'first');bdi = bdh(bid);
        bdh = bdi-bdh;
        vph = bdi-rb;
        %Prepare velocity data
        vp = v.(fn{ii});
        cp2 = find(vp.time >= start & vp.time <= stop);
        vpt = vp.time(cp2);
        x = vp.x(cp2,:);
        y = vp.y(cp2,:);
        %Rotate horiz. velocities to along & cross shore
        rot = hoffset(ii)*pi/180;
        x = x.*(ones(size(x))*cos(rot)) + ...
            y.*(ones(size(y))*sin(rot));
        y = -y.*(ones(size(y))*sin(rot)) + ...
            x.*(ones(size(x))*cos(rot));
        %Interpolate x,y (downsample) to 10Hz
        x10 = zeros(length(bdt),35);
        y10 = zeros(length(bdt),35);
        %Find unique values to avoid error with interp1
        bu = unique(bdt);
        [vu,id] = unique(vpt);
        for j = 1:35
            x10(:,j) = interp1(vu,x(id,j),bu);
            y10(:,j) = interp1(vu,y(id,j),bu);
        end
        %Prepare turbulence data
%         cp3 = find(Stat.(fn{ii}).time >= start & Stat.(fn{ii}).time <= stop);
        eps = nanmean((Stat.(fn{ii}).z1.E(bins,cp1)+Stat.(fn{ii}).z2.E(bins,cp1))./2);
        %Prepare wave energy data
        cp3 = find(wave.time2 >= start & wave.time2 <= stop);
        Ew = interp1(wave.time2(cp3),wave.E(cp3),bu);
        %Pre-process bottom distance data (see Diggle for details)
        Z = genARIMA(bh);
        %
        %
        %
        %
        
    end
end
