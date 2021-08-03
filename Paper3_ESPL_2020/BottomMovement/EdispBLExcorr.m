%Plot TKE & BLE (bed movement) for an experiment. Then, calculate
%xcorr(TKE,BLE) if sensible.
%
% Paper 3: Sediment Motion in Mangroves
%
% This is Version 1.0 of this script
%
% V 1.0: After attempting multiple methods to generate cross-correlations,
% I have landed on a method. This code produces estimates of the time-lag
% between TKE dissipation and BLE, and additionally computes a scatter plot
% of bed movement, dissipation, and bed variance.
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
load('d:\Projects\Mekong_W2015\DataAnalysis\Paper3\VTA1a_mudwvs.mat');
w1 = wave;
load('d:\Projects\Mekong_W2015\DataAnalysis\Paper3\VTA1a_fringewvs.mat');
w2 = wave;
delx = 30.61; %m
theta = 0.58; %deg
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
        E = nanmean((Stat.(fn{ii}).z1.E(bins,cp1)+Stat.(fn{ii}).z2.E(bins,cp1))./2);
        %Prepare wave energy data
        cp4 = find(w1.time2 >= start & w1.time2 <= stop);
        wef = abs((w1.F-w2.F)./(delx.*cosd(theta)));
        wef = interp1(w1.time2,wef,bu);
        %Compute var(bdh), average epsilon & bdh over interval, find
        %xcorr(epsilon,bdh) & time lags
        window = 90; %2 min averaging interval
        step = 30; %30 second step
        fs = 10;
        avt = step*fs; %samples/step
        nwin = window*fs; %samples/window
        nsamp = length(x10);
        ind = [1 avt:avt:nsamp];
        %Initialize Variables
        s_bnorm = zeros(length(ind),1);
        bvar= zeros(length(ind),1);
        bgdm = zeros(length(ind),1);
        bm = zeros(length(ind),1);
        em = zeros(length(ind),1);
        eeq = zeros(length(ind),1);
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
                %Method for computing xcorr is from Diggle:
                %1. detrend time series inputs (calc. residuals)
                eps = detrend(E(idx));
                dbi = bdh(idx);dbi(isnan(dbi)) = [];
                b = gradient(dbi);
                bgdm(j) = mean(b);
                dbi = detrend(dbi);
                bvar(j) = var(b)*10^4;
                s_bnorm(j) = bvar(j)./(U.*window); %normalized BLE
                bm(j) = bdh(idx(1));
                em(j) = mean(E(idx));
                %2. Compute xcov of residuals
                [acor,lag] = xcorr(eps,dbi,'coeff');
                %3. Compute confidence intervals:
                %3a. Compute acorr of time series
                [cor,~] = xcorr(eps,eps,'coeff');
                pxx = cor(1);
                [cor,~] = xcorr(dbi,dbi,'coeff');
                pyy = cor(1);
                %3b. Solve eqn 8.3.3
                sigs = (1+2*pxx*pyy);
                sig2 = 2*sqrt(sigs/nwin); %this is +/-95% CI
                [~,I] = max(abs(acor));
                if acor(I)>sig2 || acor(I)<-sig2
                    tD(j) = lag(I)/fs;
                else
                    tD(j) = NaN;
                end
                %Compute erosion/accretion coefficients (Eqns. 2-4, Yates
                %et al. 2009)
                a = polyfit(bdh(idx)',eps,1);
                eeq(j) = a(1)*bdh(idx(1))+a(2);
            end
        end
        %Compute erosion/accretion coefficients (Eqns. 2-4, Yates
        %et al. 2009)
        deltae = em-eeq;
        C = bgdm./(em.^(1/2).*deltae);
        Cp = C(eeq>0);Cm = C(eeq<0); %erosion/accretion coeffs
        
        
        
        tD(sT==0)=[];sT(sT==0)=[];s_bnorm(s_bnorm==0)=[];bvar(bvar==0)=[];
        figure
        sp(1) = subplot(411);
        plot(bdt,bdh,'k','linewidth',1.5)
        title('BLE')
        ylabel('m')
        sp(2) = subplot(412);
        plot(Stat.(fn{ii}).time(cp3),E,'k','linewidth',1.5)
        title('TKE dissipation rate')
        ylabel('W/kg')
        sp(3) = subplot(413);
        plot(sT,s_bnorm,'k','linewidth',1.5)
        title('Variance of BLE')
        ylabel('$\frac{\sigma_{b,int}^{2}}{U\Delta t}$',...
            'interpreter','latex')
        sp(4) = subplot(414);
        plot(sT,zeros(length(tD),1),'--k'),hold on
        plot(sT,tD,'k','linewidth',1.5)
        title('Lag time between BLE and TKE')
        ylabel('Time (s)')
        xlabel(['Time on ' datestr(bdt(1),'dd-mm-yy')])
        set([sp(1) sp(2) sp(3)],'xticklabel',[])
        set(sp,'xlim',[bdt(1) bdt(end)])
        datetickzoom('x','HH:MM','keepticks','keeplimits')
        linkaxes(sp,'x')
        
        figure
        set(gcf,'Renderer','zbuffer')
        em(bm==0)=[];bm(bm==0)=[];
        sbmin = log10(min(bvar));sbmax = log10(max(bvar));
        cax = linspace(sbmin,sbmax,length(bm));
        cl = brewermap(length(cax),'RdBu');cl = flipud(cl);
        pf = polyfit(bm,em,1);
        xs = linspace(min(bm),max(bm),length(bm));
        pv = polyval(pf,xs);
        for j = 1:length(bm)
            tmp = abs(cax-log10(bvar(j)));
            [~,I]=min(tmp);
            plot(bm(j),em(j),'o',...
                'markersize',8,...
                'markerfacecolor',cl(I,:),...
                'markeredgecolor','k'),hold on
        end
        plot(xs,pv,'-k','linewidth',1.5)
        colormap(cl)
        cb = colorbar('v');
        caxis([sbmin sbmax])
        ylabel(cb,'$log_{10}\left (\sigma_{b,int}^{2}\right )$',...
            'interpreter','latex')
        ylabel('$\left | \epsilon \right | (Wkg^{-1})$','interpreter','latex')
        xlabel('$BLE (m)$','interpreter','latex')
        set(gca,'position',[0.1 0.13 0.65 0.8],...
            'xlim',[min(bm) max(bm)])
        prettyfigures('font','Cambria','text',12,'labels',14,...
            'box',1)
        figure
        clear sp
        sp(1) = subplot(121);
        hb = herrorbar(0,mean(Cp),std(Cp),std(Cp),'s');
        set(hb,'color','k',...
            'markerfacecolor','k',...
            'markeredgecolor','k',...
            'linewidth',1.5)
        title('Erosion Coefficient')
        xlabel('m')
        sp(2) = subplot(122);
        hb = herrorbar(0,mean(Cm),std(Cm),std(Cm),'*');
        set(hb,'color','k',...
            'markerfacecolor','k',...
            'markeredgecolor','k',...
            'linewidth',1.5)
        title('Accretion Coefficient')
%         export_fig(figure(1),[sdir expname '_BLE_TKE_varBLE_lag'],'-png','-nocrop')
%         export_fig(figure(2),[sdir expname '_BLE_TKEcmap'],'-png')
%         export_fig(figure(3),[sdir expname '_Cp_Cm'],'-png')
        
    end
end




