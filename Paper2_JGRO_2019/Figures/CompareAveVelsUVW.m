%Compare mean currents of HTA1, HTA2 and VTA (VP3, VP2 and VP1). This is
%Figure 14 of Paper 2 (for now).
clear
load('D:\Projects\Mekong_W2015\DataAnalysis\Paper2\AveragedVelocities\7March2015_Vave.mat')
hta1 = Avgs.vpro3;
load('D:\Projects\Mekong_W2015\DataAnalysis\Paper2\AveragedVelocities\SingleRotation\8March2015_Vave_sngrot.mat')
hta2 = Avgs.vpro2;
load('D:\Projects\Mekong_W2015\DataAnalysis\Paper2\AveragedVelocities\SingleRotation\14March2015a_bdadj_srot.mat')
vta = Avgs.vpro1;clear Avgs
%
hta1.urms = sqrt((1/length(hta1.x))*nanmean(hta1.x(1:5,:)).^2);
hta1.vrms = sqrt((1/length(hta1.y))*nanmean(hta1.y(1:5,:)).^2);
z = nanmean((hta1.z2(1:5,:)+hta1.z2(1:5,:))./2);
hta1.wrms = sqrt((1/length(z))*z.^2);
%hta2 time has some zeros in it
ind = find(hta2.time > 1);
hta2.urms = sqrt((1/length(hta2.x(ind)))*nanmean(hta2.x(9:23,ind)).^2);
hta2.vrms = sqrt((1/length(hta2.y(ind)))*nanmean(hta2.y(9:23,ind)).^2);
z = nanmean((hta2.z2(9:23,ind)+hta2.z2(9:23,ind))./2);
hta2.wrms = sqrt((1/length(z))*z.^2);
hta2.time = hta2.time(ind);
%
vta.urms = sqrt((1/length(vta.x))*nanmean(vta.x(1:5,:)).^2);
vta.vrms = sqrt((1/length(vta.y))*nanmean(vta.y(1:5,:)).^2);
z = nanmean((vta.z2(1:5,:)+vta.z2(1:5,:))./2);
vta.wrms = sqrt((1/length(z))*z.^2);

%Plot Routine
f1 = figure(1);
set(f1,'PaperOrientation','portrait',...
    'position',[400 100   1000   400]);
set(gcf,'color','w','paperpositionmode','auto')
%
sp(1) = subplot(131);
%Horizontal
horiz = hta1.urms./hta1.vrms;
t = hta1.time;
id = find(horiz > 10 | horiz < 0.01);
horiz(id) = [];tid = setxor(1:length(t),id);
horiz = fixgaps(horiz);
hsmooth = fastsmooth(horiz,10,1,1);
plot(t(tid),horiz,...
    'Color',[0.7 0.7 0.7],...
    'Linewidth',1),hold on
plot(t(tid),hsmooth,...
    'Color','k',...
    'LineStyle',':',...
    'Linewidth',2)
%Vertical
vert = (hta1.wrms./hta1.urms)./2;
id = find(vert > 10 | vert < 0.01);
vert(id) = [];tid = setxor(1:length(t),id);
vert = fixgaps(vert);
vsmooth = fastsmooth(vert,10,1,1);
plot(t(tid),vert,...
    'Color',[0.5 0.5 0.5],...
    'Linewidth',1)
plot(t(tid),vsmooth,...
    'Color','k',...
    'LineStyle','-',...
    'Linewidth',2)
%Horizontal line at 1
plot(linspace(t(1),t(end),length(horiz)),ones(length(horiz)),...
    'Color','k','Linewidth',1)
%Plot x-scaling
tstep = datenum(0,0,0,1,0,0);
set(gca,'xlim',[t(11) t(end-9)],...
    'xtick',t(11):tstep:t(end-9),...
    'yscale','log',...
    'ylim',[10^-2 10^1],...
    'ytick',[10^-2 10^-1 10^0 10^1])
datetickzoom('x','HH:MM','keepticks','keeplimits')
%
sp(2) = subplot(132);
%Vertical
t = hta2.time;
vert = hta2.wrms./hta2.urms;
id = find(vert > 10 | vert < 0.01);
vert(id) = [];tid = setxor(1:length(t),id);
vsmooth = fastsmooth(vert,10,1,1);
plot(t(tid),vert,...
    'Color',[0.5 0.5 0.5],...
    'Linewidth',1), hold on
plot(t(tid),vsmooth,...
    'Color','k',...
    'LineStyle','-',...
    'Linewidth',2)
%Horizontal
horiz = hta2.urms./hta2.vrms;
id = find(horiz > 10 | horiz < 0.01);
horiz(id) = [];tid = setxor(1:length(t),id);
hsmooth = fastsmooth(horiz,10,1,1);
plot(t(tid),horiz,...
    'Color',[0.7 0.7 0.7],...
    'Linewidth',1)
plot(t(tid),hsmooth,...
    'Color','k',...
    'LineStyle',':',...
    'Linewidth',2)
%Horizontal line at 1
plot(linspace(t(1),t(end),length(horiz)),ones(length(horiz)),...
    'Color','k','Linewidth',1)
%Plot x-scaling
tstep = datenum(0,0,0,0,30,0);
set(gca,'xlim',[t(10) t(end-76)],...
    'xtick',t(10):tstep:t(end-76),...
    'yscale','log',...
    'ylim',[10^-2 10^1],...
    'ytick',[10^-2 10^-1 10^0 10^1])
datetickzoom('x','HH:MM','keepticks','keeplimits')
%
sp(3) = subplot(133);
%Vertical
t = vta.time;
vert = vta.wrms./vta.urms;
id = find(vert > 10 | vert < 0.01);
vert(id) = [];tid = setxor(1:length(t),id);
vert = fixgaps(vert);
vsmooth = fastsmooth(vert,10,1,1);
plot(t(tid),vert,...
    'Color',[0.5 0.5 0.5],...
    'Linewidth',1), hold on
plot(t(tid),vsmooth,...
    'Color','k',...
    'LineStyle','-',...
    'Linewidth',2)
%Horizontal
horiz = vta.urms./vta.vrms;
id = find(horiz > 10 | horiz < 0.01);
horiz(id) = [];tid = setxor(1:length(t),id);
horiz = fixgaps(horiz);
hsmooth = fastsmooth(horiz,10,1,1);
plot(t(tid),horiz,...
    'Color',[0.7 0.7 0.7],...
    'Linewidth',1)
plot(t(tid),hsmooth,...
    'Color','k',...
    'LineStyle',':',...
    'Linewidth',2)
%Horizontal line at 1
plot(linspace(t(1),t(end),length(horiz)),ones(length(horiz)),...
    'Color','k','Linewidth',1)
%Plot x-scaling
tstep = datenum(0,0,0,2,0,0);
set(gca,'xlim',[t(13) t(end)],...
    'xtick',t(13):tstep:t(end),...
    'yscale','log',...
    'ylim',[10^-2 10^1],...
    'ytick',[10^-2 10^-1 10^0 10^1])
datetickzoom('x','HH:MM','keepticks','keeplimits')

%Global Plot Settings
set(sp(1),'position',[0.08 0.15 0.26 0.7])
set(sp(2),'position',[0.4 0.15 0.26 0.7],...
    'yticklabel',[])
set(sp(3),'position',[0.72 0.15 0.26 0.7],...
    'yticklabel',[])

prettyfigures('text',13,'labels',14,'box',1)
% set(leg,'box','off')
xlabel(sp(1),['Time on ' datestr(hta1.time(1),...
    'dd-mm')])
xlabel(sp(2),['Time on ' datestr(hta2.time(1),...
    'dd-mm')])
xlabel(sp(3),['Time on ' datestr(vta.time(1),...
    'dd-mm')])
ylabel(sp(1),'u_{rms}/v_{rms},  w_{rms}/u_{rms}')
title(sp(1),'HTA1, z/h_{max} = 0.04')
title(sp(2),'HTA2, z/h_{max} = 0.33')
title(sp(3),'VTA')
figdir = 'g:\GradSchool\DataAnalysis\Paper2\WorkingFigures\Spectra\';
export_fig([figdir 'CompareRMSvels_v2'],'-pdf')
