%create an environmental conditions figure for Paper 2. Plot pressure
%signal with time divisons

clear
close all
load('D:\Projects\Mekong_W2015\Data\Vector\SW_Mudflat\V5108_030315.mat')
start = datenum(2015,03,06,13,00,00);twm = datenum(0,0,0,0,20,0);
start2 = datenum(2015,03,06,13,20,00);
stop = datenum(2015,03,06,16,40,00);
stop2 = datenum(2015,03,06,15,40,00);
ind = find(ADV.datetime >= start & ADV.datetime <= stop);
p = ADV.Pres(ind);t = ADV.datetime(ind);
pp = downsample(p,32);tt = downsample(t,32);

[~,~,~,hr,mi,~] = datevec(stop2-start2);
time = (hr*60) + mi;                                               
time = datenum(0,0,0,0,floor(time/4),0);
timee = start2:time:stop2;timee(end) = stop2;
c = brewermap(4,'YlGnBu');colormap(c)

f1 = figure(1);
set(f1,'PaperOrientation','portrait',...
    'position',[400 100   1000 500]);
set(gcf,'color','w','PaperPositionMode','auto')
%color patches
patch([start2 start2 timee(2) timee(2)],[0 2 2 0],c(1,:),...
    'EdgeColor','none');hold on
patch([timee(2) timee(2) timee(3) timee(3)],[0 2 2 0],c(2,:),...
    'EdgeColor','none');
patch([timee(3) timee(3) timee(4) timee(4)],[0 2 2 0],c(3,:),...
    'EdgeColor','none');
patch([timee(4) timee(4) timee(5) timee(5)],[0 2 2 0],c(4,:),...
    'EdgeColor','none');

%pressure signal
plot(tt,smooth(pp,1000),'LineWidth',3,'Color','k');
%divisions
line([timee(1) timee(1)],[0 2],'Color','k','LineStyle','--','LineWidth',1.5)
line([timee(2) timee(2)],[0 2],'Color','k','LineStyle',':','LineWidth',1.5)
line([timee(3) timee(3)],[0 2],'Color','k','LineStyle',':','LineWidth',1.5)
line([timee(4) timee(4)],[0 2],'Color','k','LineStyle',':','LineWidth',1.5)
line([timee(5) timee(5)],[0 2],'Color','k','LineStyle','--','LineWidth',1.5)
set(gca,...
    'XTick',start:2*twm:stop,...
    'Xlim',[start stop2+twm],...
    'Ylim',[0 1.4],...
    'Box','on',...
    'LineWidth',1.5,...
    'FontSize',14,...
    'FontName','Arial',...
    'TickDir','out')
datetickzoom('x','HH:MM:SS','keepticks','keeplimits')
xlabel('Time on 03/06/2015','FontSize',18)
ylabel('Depth (m)','FontSize',18)
fdir = 'd:\Projects\Mekong_W2015\Figures\Paper2\CurrentDirections\';
export_fig([fdir 'TideDivision'],'-pdf','-nocrop')
