%Plot Velocities, TKE as a time series and as spectra to look for
%similarities in the periodicity between wave velocity and TKE.

clear
load('D:\Projects\Mekong_W2015\DataAnalysis\Paper1\Eps_Vel_Spectra\Turbulence\HTA_07153000_30minTKE.mat')
load('D:\Projects\Mekong_W2015\DataAnalysis\Paper1\Eps_Vel_Spectra\HTA_07153000_30min.mat')
savefigdir = 'd:\Projects\Mekong_W2015\Figures\Paper1\WaveVelsTurbulence\';
plotfigs = 0;
run_ucubed = 0;
heading = 20; %degrees, angle off of transect line heading
c = [25 132 255;255 167 0]./255;
%Turbulence at VP1 should oscillate with a similar frequency as the passing
%waves because it is in front of the line of pneumatophores.

%I might need to interpolate the TKE signal to the sampling rate of the RAW
%velocities to get the correct resolution on cross-correlation or
%cospectra...
eps = (Stat.vpro1.beam1.E(15,:)+Stat.vpro1.beam2.E(15,:)+Stat.vpro1.beam3.E(15,:)+Stat.vpro1.beam4.E(15,:))./4;
eps = cmgbridge(eps,10,100,100);
time1 = Stat.vpro1.time;
x = dat.vpro1.x(:,15);y = dat.vpro1.y(:,15);
%rotate to cross-shore and along-shore
rot = heading*pi/180;
x = x.*(ones(size(x))*cos(rot)) + ...
    y.*(ones(size(y))*sin(rot));
y = -y.*(ones(size(y))*sin(rot)) + ...
    x.*(ones(size(x))*cos(rot));

time2 = dat.vpro1.time;
eps = spline(time1,eps,time2); %interpolate eps to 50Hz
eps(eps < 0) = 0; %set any negative values to zero
n = length(eps);
%for some reason the last sec of data is getting messed up by spline. Remove.
eps(n-50:n) = [];
time2(n-50:n) = [];
x(n-50:n) = [];
y(n-50:n) = [];

%03/11/16: I tried plotting time series and spectra using RAW x velocities,
%however order-of-magnitude issues may be confusing the anaylsis. Try
%calculating a*U^3, then plot this as a time series and as spectra against
%turbulence (similar order-of-magnitude)
if run_ucubed
    skipfigs = 1;
    %compute |u|^3
    ucubed = (x.^2+y.^2).^(3/2);    %'x' points 'north', 'y' points 'west'
    a = 3.5716;
    % Canopy Stats: 30x20 (WxH) box above VP1 in HTA_1
    % Max Canopy Height: 0.64
    % Delta S: 0.033387
    % Number of stems: 15
    % Mean stem diameter: 0.023878
    % a: 3.5716
    % phi: 0.078094
    x = a*ucubed;
else
    skipfigs = 0;
end

%find negative values of x & y and color them differently (log of neg is
%imaginary)
if ~skipfigs
    id = find(x<0);xneg = NaN(length(x),1);xneg(id) = x(id);
    id = find(y<0);yneg = NaN(length(x),1);yneg(id) = y(id);
    
    %raw (averaged) data plotted against one another
    start = datenum(2015,03,07,15,45,00);stop = datenum(2015,03,07,15,45,20); %snapshot of t-s
    f1 = figure(1);
    set(f1,'PaperOrientation','portrait',...
        'position',[400 200   800   500]);
    set(gcf,'color','w','PaperPositionMode','auto')
    sp(1) = subplot(211);
    p(1) = plot(time2,x,'-','color','k',...
        'LineWidth',1.5);hold on
    plot(time2,xneg,'-','color',c(1,:),...
        'linewidth',1.5)
    p(2) = plot(time2,y,'-.','color',[0.5 0.5 0.5],...
        'LineWidth',1.5);hold on
    plot(time2,yneg,'-.','color',c(2,:),...
        'linewidth',1.5)
    plot(time2,zeros(length(x),1),'color',[0.5 0.5 0.5])
    set(gca,'xlim',[start stop],'Xticklabel',[],...
        'ylim',[-0.25 0.25])
    title('Unadjusted Velocity')
    ylabel('Velocity (m/s')
    legend(p,{'x-shore';'along-shore'})
    sp(2) = subplot(212);
    plot(time2,eps,'color','k','linewidth',1.5)
    set(gca,'xlim',[start stop])
    datetickzoom('x','HH:MM:SS','keepticks','keeplimits')
    linkaxes(sp,'x')
    xlabel('Time on 07/03/2015')
    ylabel('\epsilon')
    title('Unadjusted Turbulence')
end

%Power spectra
fs = 50;
xs = detrend(x);
ys = detrend(y);
epss = detrend(eps);
step = 300; %seconds
avt = fs*step; %# samples per window
window = avt*0.2;
nfft = window*0.5; %equal to 30 second subsamples
idx = [1 avt:avt:length(xs)];
n = length(idx)-1;
% [xpsd,xfreq] = pwelch(xs,hanning(avt),avt*0.7,nfft,fs);
% [epsd,efreq] = pwelch(epss,hanning(avt),avt*0.7,nfft,fs);

xpsd = zeros(nfft/2+1,length(idx)-1);
epsd = zeros(nfft/2+1,length(idx)-1);
for i = 1:n
    xx = xs(idx(i):idx(i+1));
    epsx = epss(idx(i):idx(i+1));
    [xpsd(:,i),xfreq] = pwelch(xx,hanning(window),window*0.7,nfft,fs);
    [epsd(:,i),efreq] = pwelch(epsx,hanning(window),window*0.7,nfft,fs);
end
%magnitudes in dB:
xmag = 10*log10(xpsd);
emag = 10*log10(epsd);
tt = linspace(0,30,n);
[~,id] = max(max(xmag));min_elapsed = tt(id);

f2 = figure(2);
set(f2,'PaperOrientation','portrait',...
    'position',[400 200   800   800]);
set(gcf,'color','w','PaperPositionMode','auto')
colormap(jet)
subplot(211)
% p(1) = plot(xfreq,xmag,'-k','Linewidth',1.5);
% xlabel('Frequency (Hz)'),ylabel('dB'),title('Velocity PSD')
p(1) = surf(tt,xfreq,xmag);hold on
set(gca,'xlim',[0 max(tt)],'ylim',[0 25]),shading interp
view([90 90])
% create patch to designate the sample point for figure 3
f(1) = fill3([min_elapsed min_elapsed min_elapsed min_elapsed],[25 0 0 25],[-100 -100 max(xmag(:,id)) max(xmag(:,id))],[0.85 0.85 0.85]);
plot3(repmat(min_elapsed,length(xfreq)),xfreq,xmag(:,id),...
'-k','linewidth',1.5);
xlabel('Time (Minutes)'),ylabel('Frequency (Hz)'),zlabel('Power/Frequency (dB/Hz)')
title(['Velocity spectra, ' num2str(step) ' sec @ ' num2str(nfft/(fs)) ' sec window & 70% overlap'])
caxis([-100 -10])
cb = colorbar;
ylabel(cb,'dB')

subplot(212)
% p(2) = plot(efreq,emag,'-k','linewidth',1.5);
% xlabel('Frequency (Hz)'),ylabel('dB'),title('Turbulence PSD')
p(2) = surf(tt,efreq,emag);hold on
set(gca,'xlim',[0 max(tt)],'ylim',[0 25]),shading interp
view([90 90])
%create patch to designate the sample point for figure 3
f(2) = fill3([min_elapsed min_elapsed min_elapsed min_elapsed],[25 0 0 25],[-200 -200 max(emag(:,id)) max(emag(:,id))],[0.85 0.85 0.85]);
plot3(repmat(min_elapsed,length(efreq)),xfreq,emag(:,id),...
'-k','linewidth',1.5);
xlabel('Time (Minutes)'),ylabel('Frequency (Hz)'),zlabel('Power/Frequency (dB/Hz)')
title(['Turbulence spectra, ' num2str(step) ' sec @ ' num2str(nfft/fs) ' sec window & 70% overlap'])
set(f,'EdgeColor','none','facealpha',0.6)
caxis([-150 -50])
cb = colorbar;
ylabel(cb,'dB')

f3 = figure(3);
set(f3,'PaperOrientation','portrait',...
    'position',[400 200   800   500]);
set(gcf,'color','w','PaperPositionMode','auto')
p(1) = plot(xfreq,xmag(:,id),'r','linewidth',1.5); hold on
p(2) = plot(efreq,emag(:,id),'k','linewidth',1.5);
legend(p,{'Velocity';'Turbulence'})
xlabel('Frequency (Hz)')
ylabel('dB')
fig_time = time1(1)+datenum(0,0,0,0,min_elapsed,0);
title(['Velocity & Turbulence 60 sec power spectra at: ' datestr(fig_time,'dd-mm-yy HH:MM:SS')])

if ~skipfigs
    %try adjusting signals for lag with cross-correlation
    [acor,lag] = xcorr(xs,epss);
    [~,I] = max(abs(acor));
    lagDiff = lag(I);
    timeDiff = lagDiff/fs;
    disp(['Lag difference between velocity and turbulence: ' num2str(timeDiff) ' seconds'])
    x_adj = zeros(length(xs)+lagDiff,1);
    x_adj(1:length(xs)) = xs;
    y_adj = zeros(length(ys)+lagDiff,1);
    y_adj(1:length(ys)) = ys;
    e_adj = zeros(length(epss)+lagDiff,1);
    e_adj(lagDiff+1:end) = epss;
    e_norm = zeros(length(epss)+lagDiff,1);
    e_norm(1:length(epss)) = epss;
    t1 = (0:length(x_adj)-1)/fs;
    
    start = 900;stop = 920;
    id = find(x_adj<0);xneg = NaN(length(x_adj),1);xneg(id) = x_adj(id);
    id = find(y_adj<0);yneg = NaN(length(y_adj),1);yneg(id) = y_adj(id);
    f4 = figure(4);
    set(f4,'PaperOrientation','portrait',...
        'position',[400 200   800   500]);
    set(gcf,'color','w','PaperPositionMode','auto')
    sp(1) = subplot(211);
    p(1) = plot(t1,x_adj,'-','color','k',...
        'LineWidth',1.5);hold on
    plot(t1,xneg,'-','color',c(1,:),...
        'linewidth',1.5)
    p(2) = plot(t1,y_adj,'-.','color',[0.5 0.5 0.5],...
        'linewidth',1.5);
    plot(t1,yneg,'-.','color',c(2,:),...
        'linewidth',1.5)
    plot(t1,zeros(length(xneg),1),'Color',[0.5 0.5 0.5])
    set(gca,'xlim',[start stop],'Xticklabel',[],...
        'ylim',[-0.25 0.25])
    title('Unjusted Velocities')
    ylabel('Velocity (m/s)')
    legend(p,{'x-shore';'along-shore'})
    sp(2) = subplot(212);
    p(1) = plot(t1,e_norm,'k','linewidth',1.5); hold on
    p(2) = plot(t1,e_adj,'m','linewidth',1.5);
    legend(p,{'No Shift','Phase Shift'})
    set(gca,'xlim',[start stop])
    linkaxes(sp,'x')
    xlabel('Seconds Elapsed')
    ylabel('\epsilon')
    title('Adjusted Turbulence')
end

%CPSD
Sxy = zeros(nfft/2+1,length(idx)-1);
Cxy = zeros(nfft/2+1,length(idx)-1);
for i = 1:n
    xx = xs(idx(i):idx(i+1));
    epsx = epss(idx(i):idx(i+1));
    [Sxy(:,i),F] = cpsd(xx,epsx,hanning(window),window*0.7,nfft,fs);
    [Cxy(:,i),~] = mscohere(xx,epsx,hanning(window),window*0.7,nfft,fs);
end
f5 = figure(5);
set(f5,'PaperOrientation','portrait',...
    'position',[400 200   800   500]);
set(gcf,'color','w','PaperPositionMode','auto')
colormap(jet)
subplot(211)
p(1) = plot(F,abs(mean(Sxy,2)),'-k',...
    'LineWidth',1.5);
ylabel('Power'),xlabel('Frequency (Hz)')
title('Cross-Power Spectral Density')
%peak is at 0.2667 or 3.8 seconds. Turbulence corresponds with waves @ 3.8
%sec period? Wave peak period varies (for the time period) between 2.6 and
%4.2 seconds.
%calculate Magnitude-Squared Coherence between the two signals
subplot(212)
p(2) = plot(F,mean(Cxy,2),'-k',...
    'LineWidth',1.5);
ylabel('\gamma^2'),xlabel('Frequency (Hz)')
title('Magnitude-Squared Coherence')

if plotfigs
    prettyfigures('font','arial','text',14,'labels',16,...
        'fweight','bold','box',1','boxw',1.5)
    if ~skipfigs
        export_fig([savefigdir 'Signal'],'-png','-nocrop',f1)
        export_fig([savefigdir 'SignalShift'],'-png','-nocrop',f4)
    end
    export_fig([savefigdir 'PSD'],'-png','-nocrop',f2)
    export_fig([savefigdir 'PSDsnapshot'],'-png','-nocrop',f3)
    export_fig([savefigdir 'CPSD_MSC'],'-png','-nocrop',f5)
end
