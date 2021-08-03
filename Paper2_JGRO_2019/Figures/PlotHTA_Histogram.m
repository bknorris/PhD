%Correlates peaks in turbulence with peaks in rotated (Vectrino) velocities 
%and plots on a histogram for onshore/offshore flow. I've left the Vector
%load script in to plot wave statistics.

clear
savefigdir = 'd:\Projects\Mekong_W2015\Figures\Paper2\WaveVelsTurbulence\';
plotfigures = 0; %[off on]
bin = 5;
data = struct(); %preallocate data structure
%%%Data Pre-Processing%%%

%load vectrino/turbulence files from HTA1
datdir = 'd:\Projects\Mekong_W2015\DataAnalysis\Paper2\Eps_Vel_Spectra\';
load([datdir 'HTA_07153000_30min.mat'])
load([datdir 'Turbulence\Vertical\' 'HTA_07153000_30minTKE.mat']);
fn = fieldnames(dat);
for i = 1:length(fn) %only VP2 and VP3
    %average turbulence statistics together
    eps = mean((Stat.(fn{i}).z1.E(1:bin,:)+Stat.(fn{i}).z2.E(1:bin,:))./2);
    eps = cmgbridge(eps,10,100,100);eps(isnan(eps)) = 0;
    time1 = Stat.(fn{i}).time;time2 = dat.(fn{i}).time;
    %interpolate turbulence (@ 1 Hz) up to 50 Hz
    eps = spline(time1,eps,time2);eps(eps < 0) = 0;
    %load x,y data from dat structure, rotate to principal wave axis
    x = mean(dat.(fn{i}).x(:,1:bin),2);y = mean(dat.(fn{i}).y(:,1:bin),2);
    heading = 10; %rotate to transect; x is x-shore y is along-shore
    rot = (pi*heading)/180;
    T = [cos(rot) -sin(rot);...
        sin(rot) cos(rot)];
    vels = [x y];
    V = vels*T';
    x = V(:,1);y = V(:,2);
    data.(fn{i}).eps = eps;
    data.(fn{i}).x = detrend(x);
    data.(fn{i}).y = detrend(y);
    data.(fn{i}).time = time2;
end
clear dat Stat
%%%END Data Pre-processing%%%

%%%Find peaks in Turbulence T-S%%%
fs = 50;
pkhght = 6.5E-4; %minimum peak in TKE
pkdist = 2*fs; %peaks must be > 2sec apart
for i = 1:3
    eps = data.(fn{i}).eps;
    [pks,lc] = findpeaks(eps,...
        'minpeakheight',pkhght,...
        'minpeakdistance',pkdist);
    pkvels = data.(fn{i}).x(lc);
    data.(fn{i}).peaks = pks;
    data.(fn{i}).loc = lc;
end

%%%Plot Figures%%%
start = datenum(2015,03,07,15,47,22);
stop = datenum(2015,03,07,15,47,34);
t1 = datenum(2015,03,07,15,47,22.5);
f1 = figure(1);
set(f1,'PaperOrientation','portrait',...
    'position',[400 200   1200   500]);
set(gcf,'color','w','paperpositionmode','auto')
fn = fieldnames(data);
time = data.vpro1.time;
symb = {'o';'d';'p'};
line = {'-';'--';'-.'};
c = flipud([0.1 0.1 0.1;0.5 0.5 0.5;0.7 0.7 0.7]);
sp(1) = subplot(311);
for i = 1:3
    u = smooth(data.(fn{i}).x,10);
    lc = data.(fn{i}).loc;
    pks = data.(fn{i}).peaks;
    plot(time,u,line{i},...
        'color',c(i,:),'linewidth',1.5),hold on
    plot(time(lc),u(lc),symb{i},...
        'LineWidth',1.5,...
        'markeredgecolor','k',...
        'markerfacecolor',c(i,:),...
        'markersize',8)
end
text(t1,0.3,'Onshore'),text(t1,-0.3,'Offshore')
plot(time,zeros(length(time),1),'-k',...
    'linewidth',1)
set(gca,'xlim',[start stop],'xticklabel',[])
ylabel('Velocity (m/s)')
title('Cross-Shore Velocity, z/hc = 0.06')
sp(2) = subplot(312);
pp = zeros(3,1);
for i = 1:3
    eps = data.(fn{i}).eps;
    lc = data.(fn{i}).loc;
    pks = data.(fn{i}).peaks;
    plot(time,eps,line{i},...
        'color',c(i,:),'linewidth',1.5),hold on
    plot(time(lc),pks,symb{i},...
        'LineWidth',1.5,...
        'markeredgecolor','k',...
        'markerfacecolor',c(i,:),...
        'markersize',8)
    pp(i) = plot(time(1),0,line{i},...
        'marker',symb{i},...
        'color',c(i,:),...
        'LineWidth',1.5,...
        'markeredgecolor','k',...
        'markerfacecolor',c(i,:),...
        'markersize',8);    
end
set(gca,'xlim',[start stop])
leg = legend(pp,{'x = -10cm','x = 10cm','x = 20cm'});
title('Dissipation Rate')
yl = ylabel('$\widetilde{\epsilon} (Wkg^{-1})$');
set(yl,'Interpreter','latex','fontname','arial')
datetickzoom('x','SS','keepticks','keeplimits')
xlabel(['Time on ' datestr(time(1),'dd-mmm-yy HH:MM')])
%histogram
sp(3) = subplot(313);
ct = [1 4 7];
for i = 1:3
    cts = ct(i);
    u = data.(fn{i}).x;
    lc = data.(fn{i}).loc;
    pks = data.(fn{i}).peaks;
    pkvels = u(lc);
    s = sign(pkvels);
    neg = sum(s(:)==-1);
    pos = sum(s(:)==1);
    center = [cts-0.5 cts+0.5];
    br = bar(center,[neg pos]); hold on
    set(br,'facecolor',c(i,:),...
        'linewidth',1.5)
    text(cts-0.7,240,'Off     On')
    %UPDATE: 06/04/2017 calculate mean(eps) of peaks
    disp(fn{i})
    eps = mean(data.(fn{i}).eps(s<0));
    disp(['mean offshore eps: ' sprintf('%.2d',eps) ' W/kg'])
    eps = mean(data.(fn{i}).eps(s>0));
    disp(['mean onshore eps: ' sprintf('%.2d',eps) ' W/kg'])
end
set(gca,'ylim',[0 400],'ytick',0:100:400,...
    'xticklabel',{'','x = -10cm','','',...
    'x = 10cm','','','x = 20cm'})
ylabel('No. of Occurences')
%positioning
set(sp(1),'position',[0.08 0.56 0.42 0.3])
set(sp(2),'position',[0.08 0.15 0.42 0.3])
set(sp(3),'position',[0.575 0.1 0.4 0.83])
prettyfigures('text',13,'labels',14,'box',1)
set(leg,'position',[0.11 0.34 0.05 0.05],...
    'box','off')
savefigdir = 'd:\Projects\Mekong_W2015\Figures\Paper2\WaveVelsTurbulence\';
% export_fig([savefigdir 'Vels_eps_hist_v2'],'-pdf','-nocrop')