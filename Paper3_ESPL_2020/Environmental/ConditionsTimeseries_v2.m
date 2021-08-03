%Overview of forcing conditions fig - GM, 2018 Norris et al. 
clear, close all
dir1 = 'd:\Mekong_W2015\DataAnalysis\TOS\';
load([dir1 'Data_Fig1.mat'])
data2 = data;clear data
load([dir1 'Data_Fig1_v3.mat'])
data1 = data;
%%%Plot Routine%%%
f1 = figure(1);
set(f1,'PaperOrientation','portrait',...
    'position',[500 200 800 600],...
    'renderer','painters');
c = [207 176 126;
    60 166 74;
    4 76 41]./255;
%% Experiments
%2015
%F2F2
t3b = datenum(2015,3,5,9,0,0);t3e = datenum(2015,3,7,7,45,0);
%FSS
t4b = datenum(2015,3,9,2,0,0);t4e = datenum(2015,3,9,8,0,0);
%F2F3
t5b = datenum(2015,3,11,10,0,0);t5e = datenum(2015,3,12,12,0,0);
%% 2015
%Southwest:
%plot depth
sp(1) = subplot(411);
hold on
%Experiment Patches
min = 0;max = 3;
e(1) = patch([t3b t3e t3e t3b],[min min max max],[.9 .9 .9]);
e(2) = patch([t4b t4e t4e t4b],[min min max max],[.9 .9 .9]);
e(3) = patch([t5b t5e t5e t5b],[min min max max],[.9 .9 .9]);
set(e,'EdgeColor','none','facealpha',0.5)
%
id1 = data2.adv.WKsw.depth < 0.2;
data2.adv.WKsw.depth(id1) = NaN;
pp(1) = plot(data2.adv.WKsw.time,data2.adv.WKsw.depth,...
    'linewidth',1.5,'color',(c(1,:)));
id2 = data2.rbr.UWsw.depth < 0.2;
data2.rbr.WKsw.depth(id2) = NaN;
pp(2) = plot(data2.rbr.UWsw.time,data2.rbr.UWsw.depth,...
    'linewidth',1.5,'color',(c(2,:)));
id3 = data2.aqd.UWsw.depth < 0.26;
data2.aqd.UWsw.depth(id3) = NaN;
pp(3) = plot(data2.aqd.UWsw.time,data2.aqd.UWsw.depth,...
    'linewidth',1.5,'color',(c(3,:)));
leg = legend(pp,{'Mudflat';'Fringe';'Forest'});
text(0.02,0.8,'Depth [m]','Units','normalized')
%plot depth-av-vel
sp(2) = subplot(412);
hold on
%Experiment Patches
min = 0;max = 3;
e(1) = patch([t3b t3e t3e t3b],[min min max max],[.9 .9 .9]);
e(2) = patch([t4b t4e t4e t4b],[min min max max],[.9 .9 .9]);
e(3) = patch([t5b t5e t5e t5b],[min min max max],[.9 .9 .9]);
set(e,'EdgeColor','none','facealpha',0.5)
%
data2.adv.WKsw.vel(id1) = NaN;
plot(data2.adv.WKsw.time,data2.adv.WKsw.vel,...
    'linewidth',1.5,'color',(c(1,:))), hold on
data2.aqd.WKsw.vel(id2) = NaN;
plot(data2.aqd.WKsw.time,data2.aqd.WKsw.vel,...
    'linewidth',1.5,'color',(c(2,:)))
data2.aqd.UWsw.vel(id3) = NaN;
plot(data2.aqd.UWsw.time,data2.aqd.UWsw.vel,...
    'linewidth',1.5,'color',(c(3,:)))
text(0.02,0.8,'$|\bar{u}| \quad{[m s^{-1}]}$',...
    'interpreter','latex','Units','normalized')
%plot Hs
sp(3) = subplot(413);
hold on
%Experiment Patches
min = 0;max = 3;
e(1) = patch([t3b t3e t3e t3b],[min min max max],[.9 .9 .9]);
e(2) = patch([t4b t4e t4e t4b],[min min max max],[.9 .9 .9]);
e(3) = patch([t5b t5e t5e t5b],[min min max max],[.9 .9 .9]);
set(e,'EdgeColor','none','facealpha',0.5)
%
id = data2.adv.WKsw.Hs < 0.0025;
data2.adv.WKsw.Hs(id) = NaN;
plot(data2.adv.WKsw.time,data2.adv.WKsw.Hs,...
    'linewidth',1.5,'color',(c(1,:))), hold on
plot(data2.rbr.UWsw.time,data2.rbr.UWsw.Hs,...
    'linewidth',1.5,'color',(c(2,:)))
plot(data1.adv.WKsw.time,data1.adv.WKsw.Hs,...
    'linewidth',1.5,'color',(c(3,:)))
plot(data1.aqd.WKsw.time,data1.aqd.WKsw.Hs,...
    'linewidth',1.5,'color',(c(3,:)))
text(0.02,0.8,'H_s [m]','Units','normalized')
%plot SSC
sp(4) = subplot(414);
hold on
%Experiment Patches
min = 0;max = 3000;
e(1) = patch([t3b t3e t3e t3b],[min min max max],[.9 .9 .9]);
e(2) = patch([t4b t4e t4e t4b],[min min max max],[.9 .9 .9]);
e(3) = patch([t5b t5e t5e t5b],[min min max max],[.9 .9 .9]);
set(e,'EdgeColor','none','facealpha',0.5)
%
plot(data2.rbr.UWsw.time,data2.rbr.UWsw.SSC,...
    'linewidth',1.5,'color',(c(2,:))), hold on
id = data2.aqd.UWsw.SSC < 100;
data2.aqd.UWsw.SSC(id) = NaN;
plot(data2.aqd.UWsw.time,data2.aqd.UWsw.SSC,...
    'linewidth',1.5,'color',(c(3,:)))
text(0.02,0.8,'SSC [mg l^{-1}]','Units','normalized')
%% Global plot adjustments:
%axis lims, ticks
t1 = datenum(2015,03,03,21,30,00);
t2 = datenum(2015,03,15,06,25,22);
stp = datenum(0,0,1,0,0,0);
set(sp,'xlim',[t1 t2])
set(sp(1),'ylim',[0 2])
set(sp(2),'ylim',[0 0.75],'ytick',0:0.25:0.5)
set(sp(3),'ylim',[0 1],'ytick',0:0.33:0.66)
set(sp(4) ,'ylim',[0 3000],'ytick',0:1000:2500)
set([sp(1) sp(2) sp(3)],'xticklabel',[])

%positioning
set(sp(1),'position',[0.07 0.75 0.75 0.18])
set(sp(2),'position',[0.07 0.54 0.75 0.18])
set(sp(3),'position',[0.07 0.32 0.75 0.18])
set(sp(4),'position',[0.07 0.1 0.75 0.18])
set(leg,'position',[0.89 0.45 0.05 0.15])

%labeling
xlabel(sp(4),'Day in March, 2015')
datetick(sp(4),'x','dd','keepticks','keeplimits')
prettyfigures('text',13,'labels',14,'box',1)
sfdir = 'g:\GradSchool\DataAnalysis\Paper3\Figures\';
% export_fig([sfdir 'DataTimeSeries_v3'],'-pdf','-nocrop')
