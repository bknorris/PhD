%Plot above & intra canopy wave velocities from calculated from wavestats
%using the VPs and collocated pressure sensors. Plot HTA: VP1 wavestats
%versus the Vector, VTA: VP1 versus VP3.

clear
load('C:\Users\bkn5\Projects\Mekong_W2015\DataAnalysis\Spectra\Paper1\WaveStatsV5109_HTA.mat')
load('C:\Users\bkn5\Projects\Mekong_W2015\DataAnalysis\Paper1\HTAtimes.mat')
start = HTA.times.t1;stop = HTA.times.e1;
id = find(wvstats.time >= start & wvstats.time <= stop);
advtime = wvstats.time(id);
advwvv = wvstats.rorb(id);advwvv = runningmean(advwvv,10);
load('C:\Users\bkn5\Projects\Mekong_W2015\DataAnalysis\Spectra\Paper1\WaveStats_HTA_VTA.mat')

f2 = figure(2);
set(f2,'PaperOrientation','portrait',...
    'position',[400 200   800   400]);
set(gcf,'color','w','PaperPositionMode','auto')
sp(1) = subplot(121);
plot(advtime(15:5:end),advwvv(15:5:end),'-','Color',[0.5 0.5 0.5],'Linewidth',1.5,...
    'Marker','s','MarkerFaceColor',[1 1 1],'MarkerSize',8)
hold on
vptime = wvstats.hta.time;vpwvv = runningmean(wvstats.hta.wavev,10);
plot(vptime(35:5:end),vpwvv(35:5:end),'-','Color','k','Linewidth',1.5,...
    'Marker','^','MarkerFaceColor',[1 1 1],'MarkerSize',8)
datetick('x','HH:MM','keepticks','keeplimits')
xlabel('Time on 07/03/15','FontName','Cambria','FontSize',14)
ylabel('Wave Velocity (ms^-^1)','FontName','Cambria','FontSize',14)

sp(2) = subplot(122);
vptime = wvstats.vta.vp3.time;vpwvv = runningmean(wvstats.vta.vp3.wavev,10);
p(1) = plot(vptime(15:5:end),vpwvv(15:5:end),'-','Color',[0.5 0.5 0.5],'Linewidth',1.5,...
    'Marker','s','MarkerFaceColor',[1 1 1],'MarkerSize',8);
hold on
vptime = wvstats.vta.vp1.time;vpwvv = runningmean(wvstats.vta.vp1.wavev,10);
p(2) = plot(vptime(15:5:end),vpwvv(15:5:end),'-','Color','k','Linewidth',1.5,...
    'Marker','^','MarkerFaceColor',[1 1 1],'MarkerSize',8);
datetick('x','HH:MM','keepticks','keeplimits')
xlabel('Time on 14/03/15','FontName','Cambria','FontSize',14)

%legend
leg = legend(p,{'U_c';'U_c_,_m'});
set(leg,'FontName','Cambria','FontSize',14,'box','off',...
    'position',[0.4 0.82 0.05 0.05])

%global adjustments
da = datenum(0,0,0,0,25,0);
set(sp,'FontSize',14,'FontName','Cambria',...
    'box','on','LineWidth',1.5)
set(sp(1),'xlim',[wvstats.hta.time(1)+da wvstats.hta.time(end)],...
    'ylim',[0.1 0.35],'ytick',0.1:0.05:0.35)
set(sp(2),'xlim',[wvstats.vta.vp1.time(1) wvstats.vta.vp1.time(end)],...
    'ylim',[0.1 0.35],'ytick',0.1:0.05:0.35,'yticklabel',[])
%positioning
set(sp(1),'position',[0.1 0.14 0.4 0.82])
set(sp(2),'position',[0.55 0.14 0.4 0.82])

export_fig(['c:\Users\bkn5\Projects\Mekong_W2015\Figures\Environment\' 'CanopyWaveVels'],'-jpeg','-nocrop')