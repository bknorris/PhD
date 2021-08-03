%Plot TKE & BLE (bed movement) for an experiment. Then, calculate
%xcorr(TKE,BLE) if sensible.
%
% Paper 3: Sediment Motion in Mangroves
%
% This is Version 1.0 of this script
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear
close all
expname = 'VTA';
start = datenum('13-Mar-2015 00:00:00');
stop = datenum('15-Mar-2015 00:00:00');
wbbl = [0.0025 0.0025]; %from WaveBoundaryThickness.m
bins = 13:17;
%
dir1 = 'd:\Projects\Mekong_W2015\DataAnalysis\Paper3\BottomTrack\';
bdfiles = {'VTA1a_bdtrace.mat';'VTA2b_bdtrace.mat'};
dir2 = 'd:\Projects\Mekong_W2015\DataAnalysis\Paper3\Turbulence\';
tkfiles = {'13March2015a_VelsTKE.mat';'14March2015b_VelsTKE.mat'};
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

for i = 1:length(bdfiles)
    %Load the data
    bd = load([dir1 bdfiles{i}]);
    load([dir2 tkfiles{i}])
    fn = fieldnames(Stat);
    if strfind(expname,'VTA')
        ni = 1;
    else
        ni = 1:3;
    end
    %Initialize Figure
    for ii = ni
        f1 = figure(i);
        set(f1,'PaperOrientation','portrait',...
            'position',[400 100   800  400]);
        hold on
        cp1 = find(bd.(fn{ii}).time >= start & bd.(fn{ii}).time <= stop);
        bdt = bd.(fn{ii}).time(cp1);
        bdh = bd.(fn{ii}).bdist(cp1);
        bid = find(~isnan(bdh),1,'first');bdi = bdh(bid);
        bdh = bdi-bdh;
        vph = bdi-rb;
        %
        cp2 = find(Stat.(fn{ii}).time >= start & Stat.(fn{ii}).time <= stop);
        vpt = Stat.(fn{ii}).time(cp2);
        E = (Stat.(fn{ii}).z1.E(:,cp2)+Stat.(fn{ii}).z2.E(:,cp2))./2;
        %Interpolate E to 10 Hz
        E10 = zeros(min(size(E)),length(bdt));
        for k = 1:min(size(E))
            E10(k,:) = interp1(vpt,E(k,:),bdt);
        end
        %Locate values of E10 > bdh+wbbl
        for k = 1:length(bdh)
            id = find(vph <= bdh(k)+wbbl(i),1,'first');
            E10(id-3:end,k) = NaN;
        end
        E10 = nanmean(E10(bins,:));
        sp(1) = subplot(211);
        plot(bdt,bdh,'color',cl(ii,:),...
            'LineWidth',1.5)
        sp(2) = subplot(212);
        plot(bdt,E10,'color',cl(ii,:),...
            'LineWidth',1.5)
    end
    set(sp(1),'xticklabel',[])
    set(sp(1),'position',[0.1 0.55 0.8 0.3])
    set(sp(2),'position',[0.1 0.15 0.8 0.3])
    set(sp,'xlim',[vpt(1) vpt(end)])
    datetickzoom('x','HH:MM:SS','keepticks','keeplimits')
    xlabel(sp(2),['Time on ' datestr(vpt(1),'dd-mm-yy')])
    ylabel(sp(1),'Bed Elevation (m)'),ylabel(sp(2),'\epsilon (W/kg)')
    ttl = regexprep(bdfiles{i},'_bdtrace.mat','');
    title(sp(1),['BLE & TKE dissipation - ' ttl])
    linkaxes(sp,'x')
end
prettyfigures('text',12,'labels',14,...
    'box',1)

%Hone in on a specific time period when TKE is high; look for xcorrelation
%between the two signals
t1 = {'13-Mar-2015 16:41:31';'14-Mar-2015 13:06:52'};
t2 = {'13-Mar-2015 17:07:55';'14-Mar-2015 13:41:48'};
for i = 1:length(bdfiles)
    start = datenum(t1{i});
    stop = datenum(t2{i});
    %Load the data
    bd = load([dir1 bdfiles{i}]);
    load([dir2 tkfiles{i}])
    fn = fieldnames(Stat);
    if strfind(expname,'VTA')
        ni = 1;
    else
        ni = 1:3;
    end
    for ii = ni
        cp1 = find(bd.(fn{ii}).time >= start & bd.(fn{ii}).time <= stop);
        bdt = bd.(fn{ii}).time(cp1);
        bdh = bd.(fn{ii}).bdist(cp1);
        bid = find(~isnan(bdh),1,'first');bdi = bdh(bid);
        bdh = bdi-bdh;
        vph = bdi-rb;
        %
        cp2 = find(Stat.(fn{ii}).time >= start & Stat.(fn{ii}).time <= stop);
        vpt = Stat.(fn{ii}).time(cp2);
        E = (Stat.(fn{ii}).z1.E(:,cp2)+Stat.(fn{ii}).z2.E(:,cp2))./2;
        %Interpolate E to 10 Hz
        E10 = zeros(min(size(E)),length(bdt));
        for k = 1:min(size(E))
            E10(k,:) = interp1(vpt,E(k,:),bdt);
        end
        %Locate values of E10 > bdh+wbbl
        for k = 1:length(bdh)
            id = find(vph <= bdh(k)+wbbl(i),1,'first');
            E10(id-3:end,k) = NaN;
        end
        E10 = nanmean(E10(bins,:))';
        bdt(isnan(E10)) = [];bdh = detrend(bdh);
        bdh(isnan(E10)) = [];E10(isnan(E10)) = [];
        %XCorr
        bdh_sm = smooth(bdh,512);
        E10_sm = smooth(E10,51);
        [maxi,map] = findpeaks(bdh_sm,'minpeakdistance',600); %find maxima
        [mini,mip] = findpeaks(-bdh_sm,'minpeakdistance',600); %find minima
        pts = setxor(map,mip);
        cl = jet(length(pts));
        
        figure(1)
        for k = 1:length(pts)-1
            id = pts(k:k+1);
            [acor,lag] = xcorr(E10(id(1):id(2)),bdh_sm(id(1):id(2)));   
            [~,I] = max(abs(acor));
            lagDiff = lag(I);
            timeDiff(k) = lagDiff/10;
            
            bx = bdh(id(1):id(2));
            by = E10(id(1):id(2));
            tt = bdt(id(1):id(2));
            sp(1) = subplot(211);
            plot(tt,bx,'color',cl(k,:),...
                'LineWidth',1.5),hold on
            plot(tt(end),bx(end),'ok','markersize',10)
            sp(2) = subplot(212);
            plot(tt,by,'color',cl(k,:),...
                'LineWidth',1.5),hold on
            plot(tt(end),by(end),'ok','markersize',10)
        end
        figure
        plot(timeDiff)
        figure(1)
        linkaxes(sp,'x')
        datetickzoom('x','HH:MM:SS','keepticks','keeplimits')

