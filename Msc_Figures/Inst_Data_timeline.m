%Make a timeline for the Winter 2015 Deployment
clear
starttime = datenum(2015,03,03);
endtime = datenum(2015,03,17);

%% Dinh An waterlevel
load('c:\Users\bkn5\Projects\Mekong_W2015\DinhAnTides\dinhan.mat');

dinhanwl.wl(dinhanwl.wl == -9999) = nan;

[NAME,FREQ,TIDECON,XOUT] = t_tide(dinhanwl.wl(1:35000),'interval',.25,'start time',dinhanwl.date(1), 'latitude', 9.624);
% 1:35000 to limit to one year of data
YOUT = t_predic(dinhanwl.date,NAME,FREQ,TIDECON, 'latitude', 9.624);

p_time = starttime:datenum(0,0,0,0,30,0):endtime;
YOUT2 = t_predic(p_time,NAME,FREQ,TIDECON, 'latitude', 9.624);

%% Deployments
%F2F2
t1b = datenum(2015,3,5,9,0,0);t1e = datenum(2015,3,7,7,45,0);t1mid = ((t1b+t1e)/2)-0.2;
%FSS
t2b = datenum(2015,3,7,8,0,0);t2e = datenum(2015,3,11,9,0,0);t2mid = ((t2b+t2e)/2)-0.1;
%F2F3
t3b = datenum(2015,3,11,10,0,0);t3e = datenum(2015,3,12,12,0,0);t3mid = ((t3b+t3e)/2)-0.2;
%NE2
t4b = datenum(2015,3,13,15,30,0);t4e = datenum(2015,3,14,16,30,0);t4mid = ((t4b+t4e)/2)-0.2;

%AQDPs
%This was added 31/08/2015
%LR4
t13b = datenum(2015,3,7,12,0,0);t13e = datenum(2015,3,8,19,0,0);
%HR3
t5b = datenum(2015,3,5,12,0,0);t5e = datenum(2015,3,7,9,03,0);
t5b1 = datenum(2015,3,8,13,0,0);t5e1 = datenum(2015,3,9,16,08,0);
t5b2 = datenum(2015,3,10,13,0,0);t5e2 = datenum(2015,3,12,17,57,0);
t5b3 = datenum(2015,3,13,13,0,0);t5e3 = datenum(2015,3,14,19,40,0);
%5116
t6b = datenum(2015,3,5,12,0,0);t6e = datenum(2015,3,9,15,25,0);
t6b1 = datenum(2015,3,10,13,0,0);t6e1 = datenum(2015,3,12,20,02,0);
t6b2 = datenum(2015,3,13,13,0,0);t6e2 = datenum(2015,3,14,19,54,0);
%5117
t7b = datenum(2015,3,5,12,0,0);t7e = datenum(2015,3,9,17,58,0);
t7b1 = datenum(2015,3,10,13,0,0);t7e1 = datenum(2015,3,12,19,38,0);
t7b2 = datenum(2015,3,13,13,0,0);t7e2 = datenum(2015,3,14,19,52,0);

%Vectors
%5108
t8b = datenum(2015,3,3,18,0,0);t8e = datenum(2015,3,8,9,34,0);
t8b1 = datenum(2015,3,8,13,0,0);t8e1 = datenum(2015,3,12,10,28,0);
t8b2 = datenum(2015,3,12,13,0,0);t8e2 = datenum(2015,3,15,11,02,0);
%5109
t9b = datenum(2015,3,5,12,0,0);t9e = datenum(2015,3,9,16,27,0);
t9b1 = datenum(2015,3,10,13,0,0);t9e1 = datenum(2015,3,12,17,50,0);
t9b2 = datenum(2015,3,13,13,0,0);t9e2 = datenum(2015,3,14,19,38,0);
%VC01
t10b = datenum(2015,3,5,12,0,0);t10e = datenum(2015,3,9,16,18,0);
t10b1 = datenum(2015,3,10,13,0,0);t10e1 = datenum(2015,3,12,16,52,0);
t10b2 = datenum(2015,3,13,13,0,0);t10e2 = datenum(2015,3,14,18,46,0);
%VC101
t11b = datenum(2015,3,5,12,0,0);t11e = datenum(2015,3,9,16,20,0);
t11b1 = datenum(2015,3,10,13,0,0);t11e1 = datenum(2015,3,12,16,59,0);
t11b2 = datenum(2015,3,13,13,0,0);t11e2 = datenum(2015,3,14,18,45,0);

%Vectrinos
t12b = datenum(2015,3,5,14,0,0);t12e = datenum(2015,3,5,16,52,0);
t12b1 = datenum(2015,3,6,13,47,0);t12e1 = datenum(2015,3,6,17,01,0);
t12b2 = datenum(2015,3,7,13,51,0);t12e2 = datenum(2015,3,7,17,08,0);
t12b3 = datenum(2015,3,8,13,57,0);t12e3 = datenum(2015,3,8,19,08,0);
t12b4 = datenum(2015,3,9,2,34,0);t12e4 = datenum(2015,3,9,7,16,0);
t12b5 = datenum(2015,3,10,14,45,0);t12e5 = datenum(2015,3,10,16,39,0);
t12b6a = datenum(2015,3,11,7,31,0);t12e6a = datenum(2015,3,11,8,02,0);
t12b6b = datenum(2015,3,11,15,05,0);t12e6b = datenum(2015,3,11,16,50,0);
t12b7 = datenum(2015,3,12,7,01,0);t12e7 = datenum(2015,3,12,10,06,0);
t12b8 = datenum(2015,3,13,16,25,0);t12e8 = datenum(2015,3,13,23,55,0);
t12b9 = datenum(2015,3,14,0,03,0);t12e9 = datenum(2015,3,14,14,45,0);

%% Plot
label = {'F2F2';'FSS';'F2F3';'DPS2';'ADHR3';'AD5116';'AD5117';'V5108';'V5109';'V01';'V101';'Vectrinos';'ADLR4'};
max = 2.5;
min = -1.49;
f1 = figure(1);
set(f1,'PaperOrientation','portrait',...
    'position',[100 100   1400   800]);

ax1 = subplot(211);
%experiments
t1 = -0.2;
e(1) = patch([t1b t1e t1e t1b],[min min max max],[.9 .9 1]);t(1) = text(t1mid,t1,label{1});
e(2) = patch([t2b t2e t2e t2b],[min min max max],[.9 .9 1]);t(2) = text(t2mid,t1,label{2});
e(3) = patch([t3b t3e t3e t3b],[min min max max],[.9 .9 1]);t(3) = text(t3mid,t1,label{3});
e(4) = patch([t4b t4e t4e t4b],[min min max max],[.9 .9 1]);t(4) = text(t4mid,t1,label{4});
set(e,'EdgeColor','none')
set(t,'Color','k','FontName','Calibri','FontSize',14)

%aquadopps
a1 = 1.04; a2 = a1-0.1; a3 = a2-0.1;a4 = a1+0.1;
a(1) = line([t5b t5e],[a1 a1]);
a(2) = line([t5b1 t5e1],[a1 a1]);
a(3) = line([t5b2 t5e2],[a1 a1]);
a(4) = line([t5b3 t5e3],[a1 a1]);
a(5) = line([t6b t6e],[a2 a2]);
a(6) = line([t6b1 t6e1],[a2 a2]);
a(7) = line([t6b2 t6e2],[a2 a2]);
a(8) = line([t7b t7e],[a3 a3]);
a(9) = line([t7b1 t7e1],[a3 a3]);
a(10) = line([t7b2 t7e2],[a3 a3]);
a(11) = line([t13b t13e],[a4 a4]);
set(a,'Color',[0 0.7 1],'LineWidth',5)

%vectors
v1 = 0.7; v2 = v1-0.1; v3 = v2-0.1; v4 = v3-0.1;
v(1) = line([t8b t8e],[v1 v1]);
v(2) = line([t8b1 t8e1],[v1 v1]);
v(3) = line([t8b2 t8e2],[v1 v1]);
v(4) = line([t9b t9e],[v2 v2]);
v(5) = line([t9b1 t9e1],[v2 v2]);
v(6) = line([t9b2 t9e2],[v2 v2]);
v(7) = line([t10b t10e],[v3 v3]);
v(8) = line([t10b1 t10e1],[v3 v3]);
v(9) = line([t10b2 t10e2],[v3 v3]);
v(10) = line([t11b t11e],[v4 v4]);
v(11) = line([t11b1 t11e1],[v4 v4]);
v(12) = line([t11b2 t11e2],[v4 v4]);
set(v,'Color',[0 0.5 1],'LineWidth',5)

%vectrinos
o1 = 0.2;
o(1) = line([t12b t12e],[o1 o1]);
o(2) = line([t12b1 t12e1],[o1 o1]);
o(3) = line([t12b2 t12e2],[o1 o1]);
o(4) = line([t12b3 t12e3],[o1 o1]);
o(5) = line([t12b4 t12e4],[o1 o1]);
o(6) = line([t12b5 t12e5],[o1 o1]);
o(7) = line([t12b6a t12e6a],[o1 o1]);
o(8) = line([t12b6b t12e6b],[o1 o1]);
o(9) = line([t12b7 t12e7],[o1 o1]);
o(10) = line([t12b8 t12e8],[o1 o1]);
o(11) = line([t12b9 t12e9],[o1 o1]);
set(o,'Color',[0.5 0 0.5],'LineWidth',5)

y(1) = line([t1b t1b],[min max]);
y(2) = line([t2b t2b],[min max]);
y(3) = line([t2e t2e],[min max]);
y(4) = line([t3b t3b],[min max]);
y(5) = line([t3e t3e],[min max]);
y(6) = line([t4b t4b],[min max]);
y(7) = line([t4e t4e],[min max]);
y(8) = line([endtime endtime],[min max]);
y(9) = line([starttime endtime],[1.2 1.2]);

f(1) = text(t5b-0.82,a1,label{5});
f(2) = text(t6b-0.86,a2,label{6});
f(3) = text(t7b-0.86,a3,label{7});
f(4) = text(t8b-0.62,v1,label{8});
f(5) = text(t9b-0.76,v2,label{9});
f(6) = text(t10b-0.56,v3,label{10});
f(7) = text(t11b-0.66,v4,label{11});
f(8) = text(t12b-1.4,o1,label{12});
f(9) = text(t13b-2.8,a4,label{13});
set(y,'Color','k','LineWidth',1)
set(f,'FontSize',12,'FontName','Arial')
axis([starttime endtime -0.6 1.2])
a = gca;
set(a,...
    'XTick',[],...
    'YTick',[],...
    'YTickLabel',[],...
    'TickLength',[0 0],...
    'box','off',...
    'FontName','Calibri',...
    'FontSize',12,...
    'Linewidth',1)

ax2 = subplot(212);
e(1) = patch([t1b t1e t1e t1b],[min min max max],[.9 .9 1]);
e(2) = patch([t2b t2e t2e t2b],[min min max max],[.9 .9 1]);
e(3) = patch([t3b t3e t3e t3b],[min min max max],[.9 .9 1]);
e(4) = patch([t4b t4e t4e t4b],[min min max max],[.9 .9 1]);
set(e,'EdgeColor','none')
set(t,'Color','k','FontName','Arial','FontSize',14)
hold on
plot(p_time, YOUT2, 'linewidth', 2,'Color','k')
ylabel('\bf\itTide Level (m)','FontSize',14)
xlabel('\bf\itDate dd-mm-2015','FontSize',14)
y(1) = line([t1b t1b],[min max]);
y(2) = line([t2b t2b],[min max]);
y(3) = line([t2e t2e],[min max]);
y(4) = line([t3b t3b],[min max]);
y(5) = line([t3e t3e],[min max]);
y(6) = line([t4b t4b],[min max]);
y(7) = line([t4e t4e],[min max]);
y(8) = line([endtime endtime],[min max]);
set(y,'Color','k','LineWidth',1)
set(ax2,'YLim',[-1.5 1.5])

dis = get(ax1,'pos');
dis(2) = dis(2)-0.16;
dis(4) = dis(4)+0.1;
set(ax1,'pos',dis)
dis = get(ax2,'pos');
dis(2) = dis(2)-0.04;
dis(4) = dis(4)+0.1;
set(ax2,'pos',dis)
set(ax2,...
    'XTick',[starttime:1:endtime],...
    'TickDir','out',...
    'YTickLabel',[-1.5:0.5:1.5],...
    'box','off',...
    'FontName','Arial',...
    'FontSize',11,...
    'Linewidth',1)
datetick('x','dd/mm','keepticks','keeplimits')
set(gca,'FontSize',12,'FontName','Arial')
title = suptitle('Data Record - Mekong Winter Deployment 2015');
set(title,'FontName','Calibri','FontSize',14)
set(gcf,'color','w','PaperPositionMode','auto');
filen = 'datarecord.pdf';
filep = ['c:\Users\bkn5\Projects\Mekong_W2015\Documents\Deployment\'];
export_fig([filep filen],'-pdf','-r600')
