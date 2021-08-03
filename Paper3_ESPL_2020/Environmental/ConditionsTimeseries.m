%Overview of forcing conditions fig - GM, 2018 Norris et al. 
clear, close all
dir1 = 'd:\Projects\Mekong_W2015\DataAnalysis\TOS\';
load([dir1 'Data_Fig1.mat'])
data2 = data;clear data
load([dir1 'Data_2014.mat'])
data1 = data;
%%%Plot Routine%%%
f1 = figure(1);
set(f1,'PaperOrientation','portrait',...
    'position',[500 200 1100 500],...
    'renderer','painters');
c = [207 176 126;
    60 166 74;
    4 76 41]./255;
%% Experiments
%2014
%MA-1
t1b = datenum(2014,9,24,7,0,0);t1e = datenum(2014,9,24,16,30,0);
%F2F
t2b = datenum(2014,9,29,7,0,0);t2e = datenum(2014,9,30,10,30,0);
%2015
%F2F2
t3b = datenum(2015,3,5,9,0,0);t3e = datenum(2015,3,7,7,45,0);
%FSS
t4b = datenum(2015,3,9,2,0,0);t4e = datenum(2015,3,9,8,0,0);
%F2F3
t5b = datenum(2015,3,11,10,0,0);t5e = datenum(2015,3,12,12,0,0);
%NE2
t6b = datenum(2015,3,13,15,30,0);t6e = datenum(2015,3,14,16,30,0);

%% 2014
%Southwest
%plot depth
sp(1) = subplot(431);
hold on
%Experiment Patches
min = 0;max = 3;
e(1) = patch([t1b t1e t1e t1b],[min min max max],[.9 .9 .9]);
e(2) = patch([t2b t2e t2e t2b],[min min max max],[.9 .9 .9]);
set(e,'EdgeColor','none','facealpha',0.5)
%
plot(data1.rbr.UWsw.time,data1.rbr.UWsw.depth,...
    'linewidth',1.5,'color',(c(1,:)))
dep = data1.aqd.WKfr.depth-0.09;
id2 = dep < 0.09;dep(id2) = NaN;
plot(data1.aqd.WKfr.time,dep,...
    'linewidth',1.5,'color',(c(2,:)))
dep = data1.aqd.UWsw.depth+0.16;
id3 = dep < 0.09;dep(id3) = NaN;
plot(data1.aqd.UWsw.time,dep,...
    'linewidth',1.5,'color',(c(3,:)))
text(0.02,0.8,'Depth (m)','Units','normalized')
%plot depth-av-vel
sp(4) = subplot(434);
hold on
%Experiment Patches
min = 0;max = 3;
e(1) = patch([t1b t1e t1e t1b],[min min max max],[.9 .9 .9]);
e(2) = patch([t2b t2e t2e t2b],[min min max max],[.9 .9 .9]);
set(e,'EdgeColor','none')
%
data1.aqd.WKfr.vel(id2) = NaN;
plot(data1.aqd.WKfr.time,data1.aqd.WKfr.vel,...
    'linewidth',1.5,'color',(c(2,:))), hold on
data1.aqd.WKsw.vel(id3) = NaN;
plot(data1.aqd.UWsw.time,data1.aqd.UWsw.vel,...
    'linewidth',1.5,'color',(c(3,:)))
text(0.02,0.8,'$\langle{\bar{u}}\rangle \quad{(m s^{-1})}$',...
    'interpreter','latex','Units','normalized')
%plot Hs
sp(7) = subplot(437);
hold on
%Experiment Patches
min = 0;max = 3;
e(1) = patch([t1b t1e t1e t1b],[min min max max],[.9 .9 .9]);
e(2) = patch([t2b t2e t2e t2b],[min min max max],[.9 .9 .9]);
set(e,'EdgeColor','none')
%
data1.aqd.WKsw.Hs(id2) = NaN;
plot(data1.aqd.WKmd.time,data1.aqd.WKmd.Hs,...
    'linewidth',1.5,'color',(c(1,:))), hold on
text(0.02,0.8,'H_s (m)','Units','normalized')
plot(data1.rbr.UWsw.time,data1.rbr.UWsw.Hs,...
    'linewidth',1.5,'color',(c(2,:))), hold on
%plot SSC
sp(10) = subplot(4,3,10);
hold on
%Experiment Patches
min = 0;max = 3000;
e(1) = patch([t1b t1e t1e t1b],[min min max max],[.9 .9 .9]);
e(2) = patch([t2b t2e t2e t2b],[min min max max],[.9 .9 .9]);
set(e,'EdgeColor','none')
%
plot(data1.rbr.UWsw.time,data1.rbr.UWsw.SSC,...
    'linewidth',1.5,'color',(c(2,:))), hold on
data1.aqd.WKsw.SSC(id3) = NaN;
plot(data1.aqd.UWsw.time,data1.aqd.UWsw.SSC,...
    'linewidth',1.5,'color',(c(3,:))), hold on
text(0.02,0.8,'SSC (mg l^-^1)','Units','normalized')
%% 2015
%Southwest:
%plot depth
sp(2) = subplot(432);
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
plot(data2.adv.WKsw.time,data2.adv.WKsw.depth,...
    'linewidth',1.5,'color',(c(1,:)))
id2 = data2.rbr.UWsw.depth < 0.2;
data2.rbr.WKsw.depth(id2) = NaN;
plot(data2.rbr.UWsw.time,data2.rbr.UWsw.depth,...
    'linewidth',1.5,'color',(c(2,:)))
id3 = data2.aqd.UWsw.depth < 0.26;
data2.aqd.UWsw.depth(id3) = NaN;
plot(data2.aqd.UWsw.time,data2.aqd.UWsw.depth,...
    'linewidth',1.5,'color',(c(3,:)))
%plot depth-av-vel
sp(5) = subplot(435);
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
%plot Hs
sp(8) = subplot(438);
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
%plot SSC
sp(11) = subplot(4,3,11);
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

%Northeast:
%plot depth
sp(3) = subplot(433);
hold on
%Experiment Patches
min = 0;max = 3;
e(1) = patch([t6b t6e t6e t6b],[min min max max],[.9 .9 .9]);
set(e,'EdgeColor','none','facealpha',0.5)
%
id = data2.adv.WKne.depth < 0.13;
data2.adv.WKne.depth(id) = NaN;
pp(1) = plot(data2.adv.WKne.time,data2.adv.WKne.depth,...
    'linewidth',1.5,'color',(c(1,:)));
id = data2.rbr.UWne.depth < 0.122;
data2.rbr.UWne.depth(id) = NaN;
pp(2) = plot(data2.rbr.UWne.time,data2.rbr.UWne.depth,...
    'linewidth',1.5,'color',(c(2,:)));
id2 = data2.aqd.UWne.depth < 0.43;
data2.aqd.UWne.depth(id2) = NaN;
pp(3) = plot(data2.aqd.UWne.time,data2.aqd.UWne.depth,...
    'linewidth',1.5,'color',(c(3,:)));
leg = legend(pp,{'Mudflat';'Fringe';'Forest'});
%plot depth-av-vel
sp(6) = subplot(436);
hold on
min = 0;max = 3;
e(1) = patch([t6b t6e t6e t6b],[min min max max],[.9 .9 .9]);
set(e,'EdgeColor','none','facealpha',0.5)
%
plot(data2.adv.WKne.time,data2.adv.WKne.vel,...
    'linewidth',1.5,'color',(c(1,:)))
plot(data2.aqd.WKne.time,data2.aqd.WKne.vel,...
    'linewidth',1.5,'color',(c(2,:)))
plot(data2.aqd.UWne.time,data2.aqd.UWne.vel,...
    'linewidth',1.5,'color',(c(3,:)))
%plot Hs
sp(9) = subplot(439);
hold on
%Experiment Patches
min = 0;max = 3;
e(1) = patch([t6b t6e t6e t6b],[min min max max],[.9 .9 .9]);
set(e,'EdgeColor','none','facealpha',0.5)
%
id = data2.adv.WKne.Hs < 0.0045;
data2.adv.WKne.Hs(id) = NaN;
plot(data2.adv.WKne.time,data2.adv.WKne.Hs,...
    'linewidth',1.5,'color',(c(1,:)))
plot(data2.rbr.UWne.time,data2.rbr.UWne.Hs,...
    'linewidth',1.5,'color',(c(2,:)))
%plot SSC
sp(12) = subplot(4,3,12);
hold on
%Experiment Patches
min = 0;max = 3000;
e(1) = patch([t6b t6e t6e t6b],[min min max max],[.9 .9 .9]);
set(e,'EdgeColor','none','facealpha',0.5)
%
plot(data2.rbr.UWne.time,data2.rbr.UWne.SSC,...
    'linewidth',1.5,'color',(c(2,:)))
id = find(data2.aqd.UWne.SSC < 100);
data2.aqd.UWne.SSC(id) = NaN;
plot(data2.aqd.UWne.time,data2.aqd.UWne.SSC,...
    'linewidth',1.5,'color',(c(3,:)))

%% Global plot adjustments:
%axis lims, ticks
t1 = datenum(2014,09,23,00,00,00);
t2 = datenum(2014,10,01,00,00,00);
stp = datenum(0,0,1,0,0,0);
set([sp(1) sp(4) sp(7) sp(10)],'xtick',t1:stp:t2,'xlim',[t1 t2])
t1 = datenum(2015,03,03,21,30,00);
t2 = datenum(2015,03,15,06,25,22);
stp = datenum(0,0,1,0,0,0);
set([sp(2) sp(5) sp(8) sp(11)],'xtick',t1:stp:t2,'xlim',[t1 t2])
t1 = datenum(2015,03,13,00,00,00);
t2 = datenum(2015,03,15,06,25,22);
set([sp(3) sp(6) sp(9) sp(12)],'xtick',t1:stp:t2,'xlim',[t1 t2])

set([sp(2) sp(3) sp(5) sp(6) sp(8) sp(9) sp(11) sp(12)],...
    'yticklabel',[])
set([sp(1) sp(2) sp(3) sp(4) sp(5) sp(6) sp(7) sp(8) sp(9)],...
    'xticklabel',[])
set([sp(1) sp(2) sp(3)],'ylim',[0 2])
set([sp(4) sp(5) sp(6)],'ylim',[0 0.75],'ytick',0:0.25:0.5)
set([sp(7) sp(8) sp(9)],'ylim',[0 1],'ytick',0:0.33:0.66)
set([sp(10) sp(11) sp(12)] ,'ylim',[0 3000],'ytick',0:1000:2500)

%positioning
set(sp(1),'position',[0.06 0.75 0.29 0.17])
set(sp(2),'position',[0.38 0.75 0.34 0.17])
set(sp(3),'position',[0.74 0.75 0.13 0.17])
set(sp(4),'position',[0.06 0.54 0.29 0.17])
set(sp(5),'position',[0.38 0.54 0.34 0.17])
set(sp(6),'position',[0.74 0.54 0.13 0.17])
set(sp(7),'position',[0.06 0.33 0.29 0.17])
set(sp(8),'position',[0.38 0.33 0.34 0.17])
set(sp(9),'position',[0.74 0.33 0.13 0.17])
set(sp(10),'position',[0.06 0.12 0.29 0.17])
set(sp(11),'position',[0.38 0.12 0.34 0.17])
set(sp(12),'position',[0.74 0.12 0.13 0.17])
set(leg,'position',[0.92 0.5 0.025 0.025])

%labeling
title(sp(1),'Southwest')
title(sp(2),'Southwest')
title(sp(3),'Northeast')
xlabel(sp(10),'Day in September, 2014')
xlabel(sp(11),'Day in March, 2015')
xlabel(sp(12),'Day in March, 2015')
datetick(sp(10),'x','dd','keepticks','keeplimits')
datetick(sp(11),'x','dd','keepticks','keeplimits')
datetick(sp(12),'x','dd','keepticks','keeplimits')
prettyfigures('text',13,'labels',14,'box',1)
sfdir = 'f:\GradSchool\DataAnalysis\Paper3\Figures\';
% export_fig([sfdir 'DataTimeSeries_extra'],'-pdf','-nocrop')
