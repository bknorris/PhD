%This code plots power spectra & CPSD for all three vectrinos from the HTA
%(either day 1 or day 2). In addition, short time series of raw velocities
%and turbulence are plotted for each of the three instruments.

%%%Basic Settings
clear
close all
datdir = 'd:\Projects\Mekong_W2015\DataAnalysis\Paper1\Eps_Vel_Spectra\';
velfiles = dir([datdir '*30min.mat']);
velfiles = {velfiles.name};
tkefiles = dir([datdir 'Turbulence\' '*30minTKE.mat']);
tkefiles = {tkefiles.name};
toprocess = 1;
savefigs = 1;
timetoplot = datenum(2015,03,07,15,45,00); %time for snapshot
depname = 'HTA_1'; %for file saving conventions, also used as the fig directory name
savefigdir = ['d:\Projects\Mekong_W2015\Figures\Paper1\WaveVelsTurbulence\' depname '\'];
heading = 20; %degrees, angle off of transect line heading

%%%Load Data
for f = toprocess
    disp(['Loading ' velfiles{f}])
    disp(['Loading ' tkefiles{f}])
    load([datdir velfiles{f}])
    load([datdir 'Turbulence\' tkefiles{f}])
    
    fn = fieldnames(dat);
    c = [0 0 0;1 0 0;0 0 1];
    p3 = zeros(3,1);
    p4 = zeros(3,1);
    tsl = zeros(3,1);
    %find the shortest t-s; sometimes the instrument datasets are not
    %the same length
    for inst = 1:3
        tsl(inst) = length(dat.(fn{inst}).x);
    end
    tsl = min(tsl);
    for inst = 1:3 %loop through vectrinos
        bin = 7;
        x = dat.(fn{inst}).x(:,bin);y = dat.(fn{inst}).y(:,bin);
        
        %rotate to cross-shore and along-shore
        rot = heading*pi/180;
        x = x.*(ones(size(x))*cos(rot)) + ...
            y.*(ones(size(y))*sin(rot));
        y = -y.*(ones(size(y))*sin(rot)) + ...
            x.*(ones(size(x))*cos(rot));
        
        %average per-beam turbulence statistics into a single time-series
        eps = (Stat.(fn{inst}).beam1.E(bin,:)+Stat.(fn{inst}).beam2.E(bin,:)...
            +Stat.(fn{inst}).beam3.E(bin,:)+Stat.(fn{inst}).beam4.E(bin,:))./4;
        eps = cmgbridge(eps,10,100,100);eps(isnan(eps)) = 0;
        time1 = Stat.(fn{inst}).time;
        
        %interpolate turbulence (@ 1 Hz) up to 50 Hz
        time2 = dat.(fn{inst}).time;
        eps = spline(time1,eps,time2);
        eps(eps < 0) = 0; %set any negative values to zero
        fs = 50;
        
        %timebases are off by several milliseconds between instruments.
        %Create new timebase
%         tstep = 1/fs;
%         newtime = linspace(0,30,length(x));
        
        %crop to shortest t-s
        x = x(1:tsl);y = y(1:tsl);eps = eps(1:tsl);
        
        %detrend
        x = detrend(x);y = detrend(y);
        
        %%%Plot snapshot of velocity & turbulence & phase lag of turbulence
        %%%behind velocity
        if inst == 1
            id = find(x<0);x2 = NaN(length(x),1);x2(id) = x(id);
            id = find(y<0);y2 = NaN(length(x),1);y2(id) = y(id);
        else
            id = find(x>0);x2 = NaN(length(x),1);x2(id) = x(id);
            id = find(y>0);y2 = NaN(length(x),1);y2(id) = y(id);
        end
        [acor,lag] = xcorr(x,eps);
        [~,I] = max(abs(acor));
        lagDiff = abs(lag(I));
        timeDiff = lagDiff/fs;
        disp(['Lag difference between velocity and turbulence: ' num2str(timeDiff) ' seconds'])
        e_adj = zeros(length(eps)+lagDiff,1);
        e_adj(lagDiff+1:end) = eps;
        tid = find(time2 >= timetoplot & time2 <= timetoplot+datenum(0,0,0,0,0,30));
        
        f1 = figure;
        set(f1,'PaperOrientation','portrait',...
            'position',[400 200   800   500]);
        sp(1) = subplot(211);
        p(1) = plot(time2(tid),x(tid),'-','color','k',...
            'LineWidth',1.5);hold on
        plot(time2(tid),x2(tid),'-','color','g',...
            'linewidth',1.5)
        p(2) = plot(time2(tid),y(tid),'-.','color',[0.5 0.5 0.5],...
            'linewidth',1.5);
        plot(time2(tid),y2(tid),'-.','color','r',...
            'linewidth',1.5)
        plot(time2(tid),zeros(length(x(tid)),1),'Color',[0.5 0.5 0.5])
        set(gca,'xlim',[time2(tid(1)) time2(tid(end))],'Xticklabel',[],...
            'ylim',[-0.4 0.4],'ytick',-0.4:0.2:0.4)
        title(['Unjusted Velocities - ' upper(fn{inst}) ', ' depname])
        ylabel('Velocity (m/s)')
        legend(p,{'x-shore';'along-shore'},...
            'position',[0.87 0.83 0.05 0.05])
        sp(2) = subplot(212);
        p(1) = plot(time2(tid),eps(tid),'k','linewidth',1.5); hold on
        p(2) = plot(time2(tid),e_adj(tid),'m','linewidth',1.5);
        legend(p,{'No Shift','Phase Shift'},...
            'position',[0.87 0.36 0.05 0.05])
        set(gca,'xlim',[time2(tid(1)) time2(tid(end))])
        linkaxes(sp,'x')
        datetickzoom('x','HH:MM:SS','keepticks','keeplimits')
        xlabel(['Time on ' datestr(time2(1),'dd-mm-yyyy')])
        ylabel('\epsilon (Wkg^-^1')
        title(['Adjusted Turbulence - ' upper(fn{inst}) ', ' depname])
        
        %%%Compte power spectra of the cross-shore velocity & turbulence
        step = 300; %seconds
        avt = fs*step; %# samples per window
        window = avt*0.2;
        nfft = window*0.5; %equal to 30 second subsamples
        idx = [1 avt:avt:length(x)];
        n = length(idx)-1;
        xpsd = zeros(nfft/2+1,length(idx)-1);
        epsd = zeros(nfft/2+1,length(idx)-1);
        for i = 1:n
            xx = x(idx(i):idx(i+1));
            epsx = eps(idx(i):idx(i+1));
            [xpsd(:,i),xfreq] = pwelch(xx,hanning(window),window*0.7,nfft,fs);
            [epsd(:,i),efreq] = pwelch(epsx,hanning(window),window*0.7,nfft,fs);
        end
        %magnitudes in dB:
        xmag = 10*log10(xpsd);
        emag = 10*log10(epsd);
        tt = linspace(0,30,n);
%         [~,id] = max(max(xmag));
        id2 = 4; %lock in place so all 3 power spectra are from the same time
        min_elapsed = tt(id2);
        
        %%%Plot Figure 2: Evolutionary Power Spectra
        f2 = figure;
        set(f2,'PaperOrientation','portrait',...
            'position',[400 200   800   800]);
        set(gcf,'color','w','PaperPositionMode','auto')
        colormap(jet)
        subplot(211)
        p(1) = surf(tt,xfreq,xmag);hold on
        set(gca,'xlim',[0 max(tt)],'ylim',[0 25]),shading interp
        view([90 90])
        % create patch to designate the sample point for figure 3
        fl(1) = fill3([min_elapsed min_elapsed min_elapsed min_elapsed],[25 0 0 25],[-100 -100 max(xmag(:,id2)) max(xmag(:,id2))],[0.85 0.85 0.85]);
        plot3(repmat(min_elapsed,length(xfreq)),xfreq,xmag(:,id2),...
            '-k','linewidth',1.5);
        xlabel('Time (Minutes)'),ylabel('Frequency (Hz)'),zlabel('Power/Frequency (dB/Hz)')
        title(['Velocity spectra, ' upper(fn{inst}) ', ' depname ' ' num2str(step) ' sec @ ' num2str(nfft/(fs)) ' sec window & 70% overlap'])
        caxis([-100 -10])
        cb = colorbar;
        ylabel(cb,'dB')
        
        subplot(212)
        p(2) = surf(tt,efreq,emag);hold on
        set(gca,'xlim',[0 max(tt)],'ylim',[0 25]),shading interp
        view([90 90])
        %create patch to designate the sample point for figure 3
        fl(2) = fill3([min_elapsed min_elapsed min_elapsed min_elapsed],[25 0 0 25],[-200 -200 max(emag(:,id2)) max(emag(:,id2))],[0.85 0.85 0.85]);
        plot3(repmat(min_elapsed,length(efreq)),xfreq,emag(:,id2),...
            '-k','linewidth',1.5);
        xlabel('Time (Minutes)'),ylabel('Frequency (Hz)'),zlabel('Power/Frequency (dB/Hz)')
        title(['Turbulence spectra, ' upper(fn{inst}) ', ' depname ' ' num2str(step) ' sec @ ' num2str(nfft/fs) ' sec window & 70% overlap'])
        set(fl,'EdgeColor','none','facealpha',0.6)
        caxis([-150 -50])
        cb = colorbar;
        ylabel(cb,'dB')
        
        %%%Power Spectra from the max peak amplitude in Fig 2
        f3 = figure(3);
        set(f3,'PaperOrientation','portrait',...
            'position',[400 200   800   500]);
        set(gcf,'color','w','PaperPositionMode','auto')
        p3(inst) = plot(xfreq,xmag(:,id2),'color',c(inst,:),...
            'linewidth',1.5); hold on
        plot(efreq,emag(:,id2),'--','color',c(inst,:),...
            'linewidth',1.5);
        xlabel('Frequency (Hz)')
        ylabel('dB')
        fig_time = time1(1)+datenum(0,0,0,0,min_elapsed,0);
        title(['Velocity & Turbulence 60 sec power spectra at: ' datestr(fig_time,'dd-mm-yy HH:MM:SS')])
        
        %%%Cross-power spectra of velocity and turbulence
        Sxy = zeros(nfft/2+1,length(idx)-1);
        Cxy = zeros(nfft/2+1,length(idx)-1);
        for i = 1:n
            xx = x(idx(i):idx(i+1));
            epsx = eps(idx(i):idx(i+1));
            [Sxy(:,i),F] = cpsd(xx,epsx,hanning(window),window*0.7,nfft,fs);
            [Cxy(:,i),~] = mscohere(xx,epsx,hanning(window),window*0.7,nfft,fs);
        end
        %magnitudes in dB
        smag = 10*log10(abs(Sxy));
        f4 = figure(4);
        set(f4,'PaperOrientation','portrait',...
            'position',[400 200   800   500]);
        set(gcf,'color','w','PaperPositionMode','auto')
        colormap(jet)
        sp(1) = subplot(211);
        p4(inst) = plot(F,mean(smag,2),'color',c(inst,:),...
            'LineWidth',1.5);hold on
        ylabel('dB'),set(gca,'xticklabel',[])
        title('CPSD(Velocity, Turbulence)')
        sp(2) = subplot(212);
        plot(F,mean(Cxy,2),'color',c(inst,:),...
            'LineWidth',1.5);hold on
        ylabel('\gamma^2'),xlabel('Frequency (Hz)')
        title('Magnitude-Squared Coherence')

    end
    figure(3)
    legend(p3,fn)
    set(gca,'xlim',[0 1.5],'xscale','log')
    figure(4)
    legend(p4,fn,'position',[0.85 0.87 0.05 0.05])
    set(sp,'xlim',[0 1.5],'xscale','log')
    prettyfigures('font','arial','text',12,'labels',14,'box',1)
end
handles = findobj('type','figure');h = sort(handles);
if savefigs
    export_fig(h(1),[savefigdir 'VP1_SignalShift'],'-png')
    export_fig(h(2),[savefigdir 'VP1_PowerSpectra'],'-png')
    export_fig(h(3),[savefigdir depname '_PowerSpectra'],'-png')
    export_fig(h(4),[savefigdir depname '_CPSD'],'-png')
    export_fig(h(5),[savefigdir 'VP2_SignalShift'],'-png')
    export_fig(h(6),[savefigdir 'VP2_PowerSpectra'],'-png')
    export_fig(h(7),[savefigdir 'VP3_SignalShift'],'-png')
    export_fig(h(8),[savefigdir 'VP3_PowerSpectra'],'-png')
end
        
        

        

       