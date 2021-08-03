%Plot an 'environmental' time series overview figure for Paper 2. This
%figure includes pressure and significant wave height time series from an
%offshore and fringe velocimeter for the HTA and VTA both.

clear
load('d:\Mekong_W2015\DataAnalysis\Paper2\WaveStats_Full\V5109_070315wvs.mat')
V1 = wave;
load('d:\Mekong_W2015\DataAnalysis\Paper2\WaveStats_Full\V5109_080315wvs.mat')
V2 = wave;
load('d:\Mekong_W2015\DataAnalysis\Paper2\WaveStats_Full\V5109_100315wvs.mat')
V3 = wave;
load('d:\Mekong_W2015\DataAnalysis\Paper2\WaveStats_Full\Duet_140315wvs.mat')
V4 = wave;
load('d:\Mekong_W2015\DataAnalysis\Paper2\WaveStats_Full\V5108_030315wvs.mat')
Mud1 = wave;
load('d:\Mekong_W2015\DataAnalysis\Paper2\WaveStats_Full\V5108_080315wvs.mat')
Mud2 = wave;
load('d:\Mekong_W2015\DataAnalysis\Paper2\WaveStats_Full\V5109_130315wvs.mat')
Mud3 = wave;
%organize data files, HTA
HTA.m.time = date2doy([Mud1.time2 Mud2.time2]);
HTA.m.h = [Mud1.h-0.1 Mud2.h];
HTA.m.Hs = [Mud1.Hs Mud2.Hs];
HTA.fr.time = date2doy([V1.time2 V2.time2 V3.time2]);
HTA.fr.h = [V1.h V2.h V3.h];
HTA.fr.Hs = [V1.Hs V2.Hs V3.Hs];
%VTA
VTA.m.time = date2doy(Mud3.time2);
VTA.m.h = Mud3.h;
VTA.m.Hs = Mud3.Hs;
VTA.fr.time = date2doy(V4.time2);
VTA.fr.h = V4.h;
VTA.fr.Hs = V4.Hs;
%times
step1 = date2day(0,0,20,0,0,0);
step2 = date2day(0,0,10,0,0,0);
%experiment gaps
t1 = date2day(03,07,13,36,00,2015);e1 = date2day(03,07,17,10,00,2015);
t2 = date2day(03,08,14,15,00,2015);e2 = date2day(03,08,19,00,00,2015);
t3 = date2day(03,10,14,45,00,2015);e3 = date2day(03,10,16,40,00,2015);
t4 = date2day(03,14,06,20,00,2015);e4 = date2day(03,14,11,00,00,2015);

%plot routine
f1 = figure;
set(f1,'PaperOrientation','portrait',...
    'position',[500 100   1000   500]);
sp(1) = subplot(221);
hold on
p = patch([t1 e1 e1 t1],[0 0 2 2],[.9 .9 .9]);set(p,'EdgeColor','none')
p = patch([t2 e2 e2 t2],[0 0 2 2],[.9 .9 .9]);set(p,'EdgeColor','none')
p = patch([t3 e3 e3 t3],[0 0 2 2],[.9 .9 .9]);set(p,'EdgeColor','none')
plot(HTA.m.time,HTA.m.h,'-','LineWidth',1.5,...
    'Color','k'),hold on
plot(HTA.fr.time,HTA.fr.h,'-','LineWidth',1.5,...
    'Color',[0.5 0.5 0.5])
set(gca,'xlim',[HTA.fr.time(1) HTA.fr.time(end)],...
    'xtick',HTA.fr.time(1):step1:HTA.fr.time(end),...
    'xticklabel',[],'ylim',[0 2],'ytick',0:0.5:2)
ylabel('Tide (m)')
hold off

sp(2) = subplot(222);
hold on
p = patch([t4 e4 e4 t4],[0 0 2 2],[.9 .9 .9]);set(p,'EdgeColor','none')
plot(VTA.m.time,VTA.m.h,'-','LineWidth',1.5,...
    'Color','k'),hold on
plot(VTA.fr.time,VTA.fr.h,'-','LineWidth',1.5,...
    'Color',[0.5 0.5 0.5])
set(gca,'xlim',[VTA.fr.time(1) VTA.fr.time(end)],...
    'xtick',VTA.fr.time(1):step2:VTA.fr.time(end),...
    'xticklabel',[],'yticklabel',[],...
    'ylim',[0 2],'ytick',0:0.5:2)
hold off

sp(3) = subplot(223);
hold on
p = patch([t1 e1 e1 t1],[0 0 2 2],[.9 .9 .9]);set(p,'EdgeColor','none')
p = patch([t2 e2 e2 t2],[0 0 2 2],[.9 .9 .9]);set(p,'EdgeColor','none')
p = patch([t3 e3 e3 t3],[0 0 2 2],[.9 .9 .9]);set(p,'EdgeColor','none')
pp(1) = plot(HTA.m.time,HTA.m.Hs,'-','LineWidth',1.5,...
    'Color','k');hold on
pp(2) = plot(HTA.fr.time,HTA.fr.Hs,'-','LineWidth',1.5,...
    'Color',[0.5 0.5 0.5]);
set(gca,'xlim',[HTA.fr.time(1) HTA.fr.time(end)],...
    'xtick',HTA.fr.time(1):step1:HTA.fr.time(end),...
    'Ylim',[0 0.8])
xlabel('Day of 2015')
ylabel('H_s (m)')
hold off

sp(4) = subplot(224);
hold on
p = patch([t4 e4 e4 t4],[0 0 2 2],[.9 .9 .9]);set(p,'EdgeColor','none')
plot(VTA.m.time,VTA.m.Hs,'-','LineWidth',1.5,...
    'Color','k'),hold on
plot(VTA.fr.time,VTA.fr.Hs,'-','LineWidth',1.5,...
    'Color',[0.5 0.5 0.5])
set(gca,'xlim',[VTA.fr.time(1) VTA.fr.time(end)],...
    'xtick',VTA.fr.time(1):step2:VTA.fr.time(end),...
    'Ylim',[0 0.8],'yticklabel',[])
xlabel('Day of 2015')
hold off
leg = legend(pp,{'Mudflat';'Fringe'});
%plot adjustments
set(sp(1),'position',[0.09 0.57 0.42 0.33])
set(sp(2),'position',[0.55 0.57 0.42 0.33])
set(sp(3),'position',[0.09 0.17 0.42 0.33])
set(sp(4),'position',[0.55 0.17 0.42 0.33])
prettyfigures('text',13,'labels',14,'box',1)
set(leg,'position',[0.89 0.42 0.05 0.05],'box','off')
% export_fig('d:\Mekong_W2015\Figures\Paper2\Tide_Hs_HTAVTA','-pdf','-nocrop')




