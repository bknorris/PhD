%Plot TKE & BLE (bed movement) for an experiment. Then, calculate
%xcorr(TKE,BLE) if sensible.
%
% Paper 3: Sediment Motion in Mangroves
%
% This is Version 1.3 of this script
%
% V 1.0: xcorr of TKE and BLE was inconclusive. Sometimes TKE preceeded
% BLE, othertimes the converse.
% V 1.1: Tried calculating spectra with FFT. Required collocated pressure
% sensor to find fwc (wave-band cutoff f). New method (from Bricker &
% Montismith, 2007) uses CSD of u and w to find fwc.
% V 1.2: (See 1.1). Issue with spectral approximations of turbulence is
% determining from where to interpolate across the wave-band after
% removing wave energy in spectra. 
% V 1.3: Uses TKE dissipation calculated at 10Hz. Will run similar xcorr
% statistical associations to look for patterns in the data.
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
        for j = 1:35
            x10(:,j) = interp1(vpt,x(:,j),bdt);
            y10(:,j) = interp1(vpt,y(:,j),bdt);
        end
        %Prepare turbulence data
        cp3 = find(Stat.(fn{ii}).time >= start & Stat.(fn{ii}).time <= stop);
        E = nanmean((Stat.(fn{ii}).z1.E(bins,cp3)+Stat.(fn{ii}).z2.E(bins,cp3))./2);
        %Compute var(bdh), average epsilon & bdh over interval, find
        %xcorr(epsilon,bdh) & time lags
        window = 120; %30 sec averaging interval
        step = 60; %5 second step
        fs = 10;
        avt = step*fs; %samples/step
        nwin = window*fs; %samples/window
        nsamp = length(x10);
        ind = [1 avt:avt:nsamp];
        %Initialize Variables
        s_bnorm = zeros(length(ind),1);
        s_bint= zeros(length(ind),1);
        bm = zeros(length(ind),1);
        em = zeros(length(ind),1);
        tD = zeros(length(ind),1);
        sT = zeros(length(ind),1);
        for j = 1:length(ind)
            if abs(nsamp-ind(j)) < nwin
                continue
            else
                sT(j) = bdt(ind(j));
                idx = ind(j):ind(j)+nwin-1;
                xm = nanmean(x10(idx,bins),2);
                ym = nanmean(y10(idx,bins),2);
                U = mean(sqrt(xm.^2+ym.^2));
                eps = E(idx);
                %Calculate sigma^2_b,int (e.g. Staudt et al. 2017)
                dbi = bdh(idx);dbi(isnan(dbi)) = [];
                s_bint(j) = var(dbi);
                s_bnorm(j) = s_bint(j)./(U.*window);
                bm(j) = mean(dbi);
                em(j) = mean(eps);
                %Xcorr
                [acor,lag] = xcorr(dbi,eps);
                [~,I] = max(abs(acor));
                tD(j) = lag(I)/fs;
            end
        end
        tD(sT==0)=[];sT(sT==0)=[];s_bnorm(s_bnorm==0)=[];
        figure
        sp(1) = subplot(211);
        plot(bdt,bdh,'k','linewidth',1.5)
        title('BLE')
        ylabel('m')
        sp(2) = subplot(212);
        plot(Stat.(fn{ii}).time(cp3),E,'k','linewidth',1.5)
        title('TKE dissipation rate')
        ylabel('W/kg')
        xlabel(['Time on ' datestr(bdt(1),'dd-mm-yy')])
        set(sp(1),'xticklabel',[])
        set(sp,'xlim',[bdt(1) bdt(end)])
        datetickzoom('x','HH:MM','keepticks','keeplimits')
        figure
        sp(1) = subplot(211);
        plot(sT,s_bnorm,'k','linewidth',1.5)
        title('Variance of BLE')
        ylabel('$\frac{\sigma_{b,int}^{2}}{U\Delta t}$',...
            'interpreter','latex')
        sp(2) = subplot(212);
        plot(sT,zeros(length(tD),1),'--k'),hold on
        plot(sT,tD,'k','linewidth',1.5)
        title('Lag time between BLE and TKE')
        ylabel('Time (s)')
        xlabel(['Time on ' datestr(bdt(1),'dd-mm-yy')])
        set(sp(1),'xticklabel',[])
        set(sp,'xlim',[bdt(1) bdt(end)])
        datetickzoom('x','HH:MM','keepticks','keeplimits')
        
       
%         export_fig(figure(2),[sdir expname '_varBLExBLE_TKE'],'-png')
%         export_fig(figure(1),[sdir expname '_BLE_TKE'],'-png')

        figure
        bm(bm==0)=[];em(em==0)=[];
        cl = brewermap(length(bm),'RdBu');cl = flipud(cl);
        sbmin = log10(min(s_bnorm));sbmax = log10(max(s_bnorm));
        cax = linspace(sbmin,sbmax,length(bm));
        for j = 1:length(bm)
            tmp = abs(cax-log10(s_bnorm(j)));
            [~,I]=min(tmp);
            plot(em(j),bm(j),'o',...
                'markersize',10,...
                'markerfacecolor',cl(I,:),...
                'markeredgecolor','k'),hold on
        end
        colormap(cl)
        cb = colorbar('v');
        caxis([sbmin sbmax])
        ylabel(cb,'$log_{10}\left [\frac{\sigma_{b,int}^{2}}{U\Delta t}  \right ]$',...
            'interpreter','latex')
        xlabel('$\left | \epsilon \right | (Wkg^{-1})$','interpreter','latex')
        ylabel('$BLE (m)$','interpreter','latex')
        set(gca,'position',[0.1 0.13 0.68 0.8])
        prettyfigures('font','Cambria','text',12,'labels',14,...
            'box',1)
        export_fig(figure(1),[sdir expname '_BLE_TKE_varBLE'],'-png','-nocrop')

    end
end

            
            
            
