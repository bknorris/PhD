%plot bottom trace for analysis
clear
load('e:\Mekong_W2015\DataAnalysis\Paper3\BottomTrack\05-03-15\F2F2_1_bdtrace.mat')
f1 = figure(1);
set(f1,'PaperOrientation','portrait',...
    'position',[400 100   800   400],...
    'renderer','painters');hold on
cl = [0.7 0.7 0.7;0.4 0.4 0.4;0 0 0];
line = {'-';'--';':'};
sp(1) = plot(vpro1.time,vpro1.bdist-0.125,line{1},...
    'color',cl(1,:),...
    'linewidth',1.5);
hold on
sp(2) = plot(vpro2.time,vpro2.bdist-0.125,line{2},...
    'color',cl(2,:),...
    'linewidth',1.5);
sp(3) = plot(vpro3.time,vpro3.bdist-0.08,line{3},...
    'color',cl(3,:),...
    'linewidth',1.5);
plot(vpro3.time,zeros(1,length(vpro3.time)),'-k')
datetickzoom('x','HH:MM','keepticks','keeplimits')
set(gca,'xlim',[vpro3.time(1) vpro3.time(end)])
leg = legend(sp,{'Mudflat';'Fringe';'Forest'});
set(leg,'position',[0.18 0.35 0.05 0.05])
xlabel(['Time on ' datestr(vpro3.time(1),'dd-mm-yy')])
ylabel('Bed Level Elevation (m)')
prettyfigures('text',12,'labels',13,'box',1)
set(leg,'box','off')
sfdir = 'd:\Projects\Mekong_W2015\Figures\Paper3\Environmental\';
% export_fig([sfdir 'BottomTrace_F2F2'],'-png')