%plot wave roses for wave spreading/direction for the HTA and VTA. (Should
%be two side-by-side plots). Incl. Wave direction along the transect figure
%for Paper 2, and current direction above and inside the pneumatophores. 
%This is version 3 of this script.

clear, close all
load('d:\Projects\Mekong_W2015\DataAnalysis\Paper2\WaveStats_Full\WVStats_HTA_VTAv2.mat')
load('d:\Projects\Mekong_W2015\DataAnalysis\Paper2\WaveStats_Full\CurStats_HTA_VTAv2.mat')

%%%Plot Routine
f1 = figure(1);
set(f1,'PaperOrientation','portrait',...
    'position',[400 200   1000   700]);
set(gcf,'color','w','paperpositionmode','auto')
cmap = brewermap(35,'PuBu');
Dir1 = 360-[data.hta1.Dir; data.hta2.Dir; data.hta4.Dir];
Spd1 = sqrt(2).*[data.hta1.Uwrms data.hta2.Uwrms data.hta4.Uwrms]';
%HTA
sp(1) = subplot(321);
theta = Dir1*pi/180;
for i = 1:35
    th = theta(:,i);
    pr = rose(th,50);if i == 1, hold on, end
    hpatch = patch(get(pr,'XData'),get(pr,'YData'),cmap(i,:),'EdgeColor','none');
    set(pr,'color',cmap(i,:))
    hold on
end
caxis(1./[1.1 0.1])
set(gca,'View',[-90 90],'YDir','reverse');

Mud = 360-mean([data.sw1.Dir(:,1:5); data.sw2.Dir(:,1:5)],2);
Mud(Mud == 360) = NaN;
Mud = downsample(Mud,2);
Mudt = date2doy([data.sw1.time2 data.sw2.time2])';
Mudt = downsample(Mudt,2);
Fr = mean(Dir1(:,1:5),2);
Fr(Fr == 360) = NaN;Fr(Fr < 330) = NaN;
Fr = downsample(Fr,2);
Frt = date2doy([data.hta1.time2 data.hta2.time2 data.hta4.time2])';
Frt = downsample(Frt,2);

sp(2) = subplot(323);
tx = 66:70;
plot(tx,350*ones(length(tx),1),...
    'linewidth',1.5,...
    'Color','k'), hold on
plot(Mudt,Mud,'o',...
    'linewidth',1,...
    'markerfacecolor',[0.3 0.3 0.3],...
    'color',[0.1 0.1 0.1],...
    'markersize',5)
plot(Frt,Fr,'^',...
    'linewidth',1,...
    'markerfacecolor',[0.7 0.7 0.7],...
    'color',[0.6 0.6 0.6],...
    'markersize',5)
set(gca,'xlim',[Mudt(1)-0.2 Mudt(end)+0.2],...
    'ylim',[270 360],'ytick',270:30:360)

Up = [sw.day1.udir sw.day2.udir sw.day3.udir];
Low = [sw.day1.ldir sw.day2.ldir sw.day3.ldir];
tt = date2doy([sw.day1.time sw.day2.time sw.day3.time])';
sp(3) = subplot(325);
plot(tx,350*ones(length(tx),1),...
    'linewidth',1.5,...
    'Color','k'), hold on
plot(tt,Up,'o',...
    'linewidth',1,...
    'markerfacecolor',[0.3 0.3 0.3],...
    'color',[0.1 0.1 0.1],...
    'markersize',5)
% plot(tt,Low,'^',...
%     'linewidth',1,...
%     'markerfacecolor',[0.7 0.7 0.7],...
%     'color',[0.6 0.6 0.6],...
%     'markersize',5)
set(gca,'xlim',[tt(1)-0.2 tt(end)+0.2],...
    'ylim',[270 360],'ytick',270:30:360)

ylabel('Direction (deg)')
xlabel('Day of 2015')
break
%VTA
sp(3) = subplot(222);
Dir2 = 360-data.vta.Dir;
theta = Dir2*pi/180;
cmap = brewermap(11,'PuBu');
for i = 1:11
    th = theta(:,i);
    pr = rose(th,50);if i == 1, hold on, end
    hpatch = patch(get(pr,'XData'),get(pr,'YData'),cmap(i,:),'EdgeColor','none');
    set(pr,'color',cmap(i,:))
    hold on
end
caxis(1./[1.1 0.1])
set(gca,'View',[-90 90],'YDir','reverse');
cmap = brewermap([],'PuBu');
colormap(cmap),cb = colorbar;
ylabel(cb,'Wave Period (s)')

Mud = 360-mean(data.ne1.Dir(:,1:5),2);
Mud(Mud == 360) = NaN;
Mudt = date2doy(data.ne1.time2)';
Fr = mean(Dir2(:,1:5),2);
Fr(Fr == 360) = NaN;
Frt = date2doy(data.vta.time2)';

sp(4) = subplot(224);
plot(Mudt,284*ones(length(Mudt),1),...
    'linewidth',1.5,...
    'Color','k'), hold on
pp(1) = plot(Mudt,Mud,'o',...
    'linewidth',1,...
    'markerfacecolor',[0.3 0.3 0.3],...
    'color',[0.1 0.1 0.1],...
    'markersize',5);
pp(2) = plot(Frt,Fr,'^',...
    'linewidth',1,...
    'markerfacecolor',[0.7 0.7 0.7],...
    'color',[0.6 0.6 0.6],...
    'markersize',5);
set(gca,'xlim',[Mudt(1) Mudt(end)],...
    'ylim',[270 360],'ytick',270:30:360)
xlabel('Day of 2015')
leg = legend(pp,{'Mudflat';'Fringe'});
%%%Plot Adjustments%%%
set(sp(1),'position',[0.1 0.52 0.4 0.4])
set(sp(3),'position',[0.45 0.52 0.4 0.4])
set(sp(2),'position',[0.1 0.15 0.35 0.25])
set(sp(4),'position',[0.52 0.15 0.35 0.25])
set(cb,'position',[0.815 0.52 0.015 0.4])
prettyfigures('text',13,'labels',14,'box',1)
set(leg,'box','off')
% export_fig('d:\Mekong_W2015\Figures\Paper2\WaveDir_HTA_VTA_v2','-pdf','-nocrop')

