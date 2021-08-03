clear, close all
load('d:\Projects\Mekong_W2015\DataAnalysis\Paper2\WaveStats_Full\UcUwUpLow_HTA_VTA')
figdir = 'f:\GradSchool\DataAnalysis\Paper2\WorkingFigures\Environmental\';

%%%Plot Routine
f1 = figure(1);
set(f1,'PaperOrientation','portrait',...
    'position',[400 200   1000   500]);
set(gcf,'color','w','paperpositionmode','auto')

%% VTA
xlim = [datenum(2015,03,14,04,00,00) datenum(2015,03,14,12,00,00)];
ylim = [-0.25 0.25];
%Uw
sp(3) = subplot(222);hold on
set(sp(3),'xticklabel',[],...
    'ytick',[-0.2 0 0.2],...
    'position',[0.7 0.58 0.25 0.35])
x = ne.day4.utime;
mag = ne.day4.uwspd;
cmap = brewermap(round(1.5*length(mag)),'PuBu');
dir = ne.day4.uwdir;
h = sticks(x,mag,dir,xlim,ylim,0.25);
for i = 1:length(x)
set(h(i),'linewidth',1.5,'color',cmap(i,:)),end
set(h(end),'linewidth',1.5,'color','k')
text(0.1,0.2,'0.25 ms^{-1}','units','normalized')
plot(linspace(xlim(1),xlim(2),10),zeros(1,10),'-k')

%Uc
sp(4) = subplot(224);hold on
set(sp(4),...
    'ytick',[-0.1 0 0.1],...
    'position',[0.7 0.15 0.25 0.35])
mag = ne.day4.ucspd;
dir = ne.day4.ucdir;
ylim = [-0.15 0.15];
h = sticks(x,mag,dir,xlim,ylim,0.1);
for i = 1:length(x)
set(h(i),'linewidth',1.5,'color',cmap(i,:)),end
set(h(end),'linewidth',1.5,'color','k')
text(0.1,0.2,'0.1 ms^{-1}','units','normalized')
plot(linspace(xlim(1),xlim(2),10),zeros(1,10),'-k')
set([sp(3) sp(4)],'xtick',xlim(1):datenum(0,0,0,4,0,0):xlim(2))
datetick('x','dd HH:MM','keepticks','keeplimits')
%% HTA
xlim = [datenum(2015,03,07,12,00,00) datenum(2015,03,11,00,00,00)];
%Uw
sp(1) = subplot(221);hold on
set(sp(1),'xticklabel',[],...
    'ytick',[-1.5 0 1.5],...
    'position',[0.12 0.58 0.52 0.35])
ylim = [-2 2];
x = [sw.day1.utime sw.day2.utime sw.day3.utime];
mag = [sw.day1.uwspd sw.day2.uwspd sw.day3.uwspd.*1.5];
[~,f,l] = cmgidgaps(mag);
cmap = [zeros(length(f(1):l(1)),3); brewermap(length(l(1)+1:f(2)-1),'PuBu');...
    zeros(length(f(2):l(2)),3); brewermap(length(l(2)+1:f(3)-1),'PuBu');...
    zeros(length(f(3):l(3)),3); brewermap(length(l(3)+1:f(4)-1),'PuBu');...
    zeros(length(f(4):l(4)),3); brewermap(length(l(4)+1:f(5)-1),'PuBu');...
    zeros(length(f(5):l(5)),3); brewermap(length(l(5)+1:f(6)-1),'PuBu');...
    zeros(length(f(6):l(6)),3); brewermap(length(l(6)+1:f(7)-1),'PuBu');...
    zeros(length(f(7):l(7)),3)];
dir = [sw.day1.uwdir sw.day2.uwdir sw.day3.uwdir];
h = sticks(x,mag,dir,xlim,ylim,1);
for i = 1:length(x)
set(h(i),'linewidth',1.5,'color',cmap(i,:)),end
set(h(end),'linewidth',1.5,'color','k')
text(0.1,0.2,'1 ms^{-1}','units','normalized')
plot(linspace(xlim(1),xlim(2),10),zeros(1,10),'-k')
%Uc
sp(2) = subplot(223);hold on
set(sp(2),...
    'ytick',[-0.4 0 0.4],...
    'position',[0.12 0.15 0.52 0.35])
ylim = [-0.5 0.5];
mag = [sw.day1.ucspd sw.day2.ucspd sw.day3.ucspd*1.5];
dir = [sw.day1.ucdir sw.day2.ucdir sw.day3.ucdir];
h = sticks(x,mag,dir,xlim,ylim,0.25);
for i = 1:length(x)
set(h(i),'linewidth',1.5,'color',cmap(i,:)),end
set(h(end),'linewidth',1.5,'color','k')
text(0.1,0.2,'0.25 ms^{-1}','units','normalized')
plot(linspace(xlim(1),xlim(2),10),zeros(1,10),'-k')
%Plot Adjustments
set([sp(1) sp(2)],'xtick',xlim(1):datenum(0,0,1,0,0,0):xlim(2))
datetick('x','dd HH:MM','keepticks','keeplimits')
%% Final Adjustments
title(sp(1),'HTA')
title(sp(3),'VTA')
ylabel(sp(1),'u_w (ms^{-1})')
ylabel(sp(2),'u_c (ms^{-1})')
xlabel(sp(2),'Day in March, 2015')
xlabel(sp(4),'Day in March, 2015')
prettyfigures('text',12,'labels',14,'box',1)
export_fig([figdir 'UwUcWhisks_v2'],'-pdf','-nocrop')
