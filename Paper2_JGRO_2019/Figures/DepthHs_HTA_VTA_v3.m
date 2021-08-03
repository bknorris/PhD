%Plot an 'environmental' time series overview figure for Paper 2. This
%figure includes pressure and significant wave height time series from an
%offshore and fringe velocimeter for the HTA and VTA both.
%This is version 2 of this script. Edits: added third row for peak period
%(Tr) at Julia's request, 16/09/2017 

clear, close all
figdir = 'f:\GradSchool\DataAnalysis\Paper2\WorkingFigures\Environmental\';
load('d:\Projects\Mekong_W2015\DataAnalysis\Paper2\WaveStats_Full\V5109_070315wvs.mat')
V1 = wave;
load('d:\Projects\Mekong_W2015\DataAnalysis\Paper2\WaveStats_Full\V5109_080315wvs.mat')
V2 = wave;
load('d:\Projects\Mekong_W2015\DataAnalysis\Paper2\WaveStats_Full\V5109_100315wvs.mat')
V3 = wave;
load('d:\Projects\Mekong_W2015\DataAnalysis\Paper2\WaveStats_Full\AD5116_15March2015wvs.mat')
V4 = wave;
%organize data files, HTA
HTA.fr.time = [V1.time2 V2.time2 V3.time2];
HTA.fr.h = [V1.h V2.h V3.h];
HTA.fr.tr = [V1.Tr V2.Tr V3.Tr];
HTA.fr.Hs = [V1.Hs V2.Hs V3.Hs];
%VTA
VTA.fr.time = V4.time2;
VTA.fr.h = V4.h;
VTA.fr.tr = V4.Tr; %odd behavior in Tr, smooth for display
VTA.fr.Hs = V4.Hs;
%times
step1 = date2day(0,0,20,0,0,0);
step2 = date2day(0,0,10,0,0,0);
%experiment gaps
t1 = datenum(2015,03,07,13,36,00);e1 = datenum(2015,03,07,17,10,00);
t2 = datenum(2015,03,08,14,15,00);e2 = datenum(2015,03,08,19,00,00);
t3 = datenum(2015,03,10,14,45,00);e3 = datenum(2015,03,10,16,40,00);
t4 = datenum(2015,03,14,04,40,00);e4 = datenum(2015,03,14,10,40,00);

%plot routine
f1 = figure;
set(f1,'PaperOrientation','portrait',...
    'position',[500 100   1000   500]);
sp(1) = subplot(321);
hold on
p = patch([t1 e1 e1 t1],[0 0 2 2],[.9 .9 .9]);set(p,'EdgeColor','none')
p = patch([t2 e2 e2 t2],[0 0 2 2],[.9 .9 .9]);set(p,'EdgeColor','none')
p = patch([t3 e3 e3 t3],[0 0 2 2],[.9 .9 .9]);set(p,'EdgeColor','none')
plot(HTA.fr.time,HTA.fr.h,'-','LineWidth',1.5,...
    'Color','k')
set(gca,'xlim',[HTA.fr.time(1) HTA.fr.time(end)],...
    'xtick',HTA.fr.time(1):step1:HTA.fr.time(end),...
    'xticklabel',[],'ylim',[0 2],'ytick',0:0.5:2)
ylabel('Water Depth (m)')
hold off

sp(2) = subplot(322);
hold on
p = patch([t4 e4 e4 t4],[0 0 2 2],[.9 .9 .9]);set(p,'EdgeColor','none')
plot(VTA.fr.time,VTA.fr.h,'-','LineWidth',1.5,...
    'Color','k')
set(gca,'xlim',[VTA.fr.time(1) VTA.fr.time(end)],...
    'xtick',VTA.fr.time(1):step2:VTA.fr.time(end),...
    'xticklabel',[],'yticklabel',[],...
    'ylim',[0 2],'ytick',0:0.5:2)
hold off

sp(3) = subplot(323);
hold on
p = patch([t1 e1 e1 t1],[0 0 6 6],[.9 .9 .9]);set(p,'EdgeColor','none')
p = patch([t2 e2 e2 t2],[0 0 6 6],[.9 .9 .9]);set(p,'EdgeColor','none')
p = patch([t3 e3 e3 t3],[0 0 6 6],[.9 .9 .9]);set(p,'EdgeColor','none')
pp(2) = plot(HTA.fr.time,HTA.fr.tr,'-','LineWidth',1.5,...
    'Color','k');
set(gca,'xlim',[HTA.fr.time(1) HTA.fr.time(end)],...
    'xtick',HTA.fr.time(1):step1:HTA.fr.time(end),...
    'Ylim',[1 5],'ytick',1:2:5,'xticklabel',[])
ylabel('Peak Period (s)')
hold off

sp(4) = subplot(324);
hold on
p = patch([t4 e4 e4 t4],[0 0 6 6],[.9 .9 .9]);set(p,'EdgeColor','none')
pp(2) = plot(VTA.fr.time,VTA.fr.tr,'-','LineWidth',1.5,...
    'Color','k');
set(gca,'xlim',[VTA.fr.time(1) VTA.fr.time(end)],...
    'xtick',VTA.fr.time(1):step2:VTA.fr.time(end),...
    'Ylim',[1 5],'ytick',1:2:5,'xticklabel',[],'yticklabel',[])
hold off

sp(5) = subplot(325);
hold on
p = patch([t1 e1 e1 t1],[0 0 2 2],[.9 .9 .9]);set(p,'EdgeColor','none')
p = patch([t2 e2 e2 t2],[0 0 2 2],[.9 .9 .9]);set(p,'EdgeColor','none')
p = patch([t3 e3 e3 t3],[0 0 2 2],[.9 .9 .9]);set(p,'EdgeColor','none')
pp(2) = plot(HTA.fr.time,HTA.fr.Hs,'-','LineWidth',1.5,...
    'Color','k');
set(gca,'xlim',[HTA.fr.time(1) HTA.fr.time(end)],...
    'xtick',HTA.fr.time(1):step1:HTA.fr.time(end),...
    'Ylim',[0 0.8])
xlabel('Day in March, 2015')
ylabel('H_s (m)')
hold off

sp(6) = subplot(326);
hold on
p = patch([t4 e4 e4 t4],[0 0 2 2],[.9 .9 .9]);set(p,'EdgeColor','none')
plot(VTA.fr.time,VTA.fr.Hs,'-','LineWidth',1.5,...
    'Color','k')
set(gca,'xlim',[VTA.fr.time(1) VTA.fr.time(end)],...
    'xtick',VTA.fr.time(1):step2:VTA.fr.time(end),...
    'Ylim',[0 0.8],'yticklabel',[])
xlabel('Day in March, 2015')
hold off

%plot adjustments
set(sp(1),'position',[0.09 0.7 0.42 0.22])
set(sp(2),'position',[0.55 0.7 0.4 0.22])
set(sp(3),'position',[0.09 0.42 0.42 0.22])
set(sp(4),'position',[0.55 0.42 0.4 0.22])
set(sp(5),'position',[0.09 0.14 0.42 0.22])
set(sp(6),'position',[0.55 0.14 0.4 0.22])
xlim1 = [datenum(2015,03,07,12,00,00) datenum(2015,03,11,00,00,00)];
xlim2 = [datenum(2015,03,13,14,00,00) datenum(2015,03,14,14,00,00)];
set([sp(1) sp(3) sp(5)],...
    'xlim',[xlim1(1) xlim1(2)],...
    'xtick',xlim1(1):datenum(0,0,1,0,0,0):xlim1(2))
set([sp(2) sp(4) sp(6)],...
    'xlim',[xlim2(1) xlim2(2)],...
    'xtick',xlim2(1):datenum(0,0,0,8,0,0):xlim2(2))
datetick(sp(5),'x','dd HH:MM','keepticks','keeplimits')
datetick(sp(6),'x','dd HH:MM','keepticks','keeplimits')
prettyfigures('text',12,'labels',14,'box',1)
% export_fig([figdir 'HTrHsTimeseries'],'-pdf','-nocrop')




