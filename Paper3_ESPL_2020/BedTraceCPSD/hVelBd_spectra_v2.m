%Plot TKE & BLE (bed movement) for an experiment. Then, calculate
%CPSD(h and x and y,BLE)
%
% Paper 3: Sediment Motion in Mangroves
%
% This is Version 2.0 of this script
%
% V 2.0: 
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
%Pressure sensor
load('D:\Projects\Mekong_W2015\Data\RBR\Duet\DPS2\Duet_140315.mat')
ll = sin(RBR.Metadata.latitude/57.29578)^2;
zp = str2double(regexp(RBR.Metadata.hab,'^\S+','match'))/1000;
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
t1 = {'13-Mar-2015 16:41:31'};
t2 = {'13-Mar-2015 17:07:55'};
for i = 1
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
        for j = 1:35
            x10(:,j) = interp1(vpt,x(:,j),bdt);
            y10(:,j) = interp1(vpt,y(:,j),bdt);
        end
        %Prepare pressure data, downsample to 10Hz
        cp3 = find(RBR.Datetime >= start & RBR.Datetime <= stop);
        pt = RBR.Datetime(cp3);
        pres = RBR.SeaPres(cp3);
        p10 = interp1(pt,pres,bdt);
        %Compute depth
        g = zeros(length(p10),1);h = zeros(length(p10),1);
        for j = 1:length(p10)
            g(j,:) = 9.780318*(1.0+(5.2788E-3+2.36E-5*ll)*ll)+1.092E-6*p10(j,:);
            h(j,:) = ((((-1.82E-15*p10(j,:)+2.279E-10)*p10(j,:)-2.2512E-5)*p10(j,:)+9.72659)*p10(j,:))/g(j,:);
        end
        
        %Set up spectra settings
        fs = 10;
        win = 600; %seconds
        step = 300; %seconds
        %load data to run processing on into structure
        data.h = h;
        data.bd = bdh;
        data.x = x10;
        data.y = y10;
        w = 600;
        %run spectra
        vars = compute_spectra(bdt,data,win,step,w,fs);
        %plot
        fname = [expname '_' fn{ii}];
        plot_spec_results(fname,vars)
    end
end
