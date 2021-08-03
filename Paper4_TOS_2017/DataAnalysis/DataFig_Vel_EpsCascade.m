%Data figure for TOS: Velocity comparison between Control and DPS2, plus
%epsilon & SSC cascade
clear
%Load Steve's data from Control & DPS2
load('D:\Projects\Mekong_W2015\DataAnalysis\ToS\dat4julia.mat')
fr = vel4julia;
time1 = t4julia;
to(1) = t0_4julia;
load('D:\Projects\Mekong_W2015\DataAnalysis\ToS\dat4julia_flats.mat')
mud = vel4julia;
time2 = t4julia;
to(2) = t0_4julia;
clear t0_4julia vel4julia t4julia
%%%
tod = doy2date(to,[2015 2014]);
load('d:\Projects\Mekong_W2015\DataAnalysis\Paper2\WaveStats_Full\Duet_140315wvs')
wvstats = wave;
load('D:\Projects\Mekong_F2014\DataAnalysis\HR3wvs.mat')

%find depth & Hs at the time of Steve's samples for CNTRL & DPS2
tmp = abs(wvstats.time2-tod(1));
[~,id] = min(tmp);
depth(1) = wvstats.h(id);
Hs(1) = wvstats.Hs(id);
tmp = abs(wave.time-tod(2));
[~,id] = min(tmp);
depth(2) = wave.p(id);
tmp = abs(wave.time2-tod(2));
[~,id] = min(tmp);
Hs(2) = wave.Hs(id);
%%%
disp(['DPS2 Hs: ' num2str(Hs(1)) ' m'])
disp(['CNTRL Hs: ' num2str(Hs(2)) ' m'])
disp(['DPS2 depth: ' num2str(depth(1)) ' m'])
disp(['CNTRL depth: ' num2str(depth(2)) ' m'])

%Load Epsilon data from DPS2, VP1
load('D:\Projects\Mekong_W2015\DataAnalysis\Paper1\TKE\Vertical\HTA_2TKE.mat');
start(1) = datenum(2015,03,08,14,15,00);stop(1) = datenum(2015,03,08,19,15,00);
start(2) = datenum(2015,03,08,15,00,00);stop(2) = datenum(2015,03,08,17,00,00);
start(3) = datenum(2015,03,08,16,05,00);stop(3) = datenum(2015,03,08,16,15,00);
start(4) = datenum(2015,03,08,16,06,00);stop(4) = datenum(2015,03,08,16,11,00);
eps = log10((Stat.vpro1.z1.E+Stat.vpro1.z2.E)./2);
rb = 0.07-0.04-linspace(0,0.03,30);Stime = Stat.vpro1.time;
b = zeros(4,1);e = zeros(4,1);
for i = 1:4
    b(i) = find(Stime >= start(i) & Stime <= stop(i),1,'first');
    e(i) = find(Stime >= start(i) & Stime <= stop(i),1,'last');
end

%Load SSC data from Aaron, near HTA2
load('D:\Projects\Mekong_W2015\DataAnalysis\TOS\RBR_SW_2015_HR.mat')
%pad burts with NaNs
nans = NaN(1552,1792);padturb = [RBR_Data.turb_burst; nans];
turb = reshape(padturb,(3600*1792),1);
SSC = 1.6252.*turb-68.382;
rbrtime = linspace(RBR_Data.time(1),RBR_Data.time(end),length(turb));
id = find(rbrtime >= start(1) & rbrtime <= stop(1));
SSC = SSC(id);
rbrtime = rbrtime(id);
%average into 1 minute blocks
avt = 6*0.5; %1 sec step
nwin = 6*20; %20 sec window
nsamp = length(SSC);
ind = [1 avt:avt:nsamp];
SSC2 = zeros(length(ind),1);
time2 = zeros(length(ind),1);
for i = 1:length(ind)
    if abs(nsamp-ind(i)) < nwin  %skip the last few indexes approaching the end of the t-s
        continue
    else
        idx = ind(i):ind(i)+nwin-1;
    end
    time2(i) = rbrtime(idx(1));
    SSC2(i) = nanmean(SSC(idx));
end
time2(time2 == 0) = [];
SSC2(SSC2 == 0) = [];

r = zeros(4,1);q = zeros(4,1);
for i = 1:4
    r(i) = find(time2 >= start(i) & time2 <= stop(i),1,'first');
    q(i) = find(time2 >= start(i) & time2 <= stop(i),1,'last');
end
%Plot Routine
%Will split upper half and lower half into separate figures and combine
%them in Illustrator
% f1 = figure(1);
% set(f1,'PaperOrientation','portrait',...
%     'position',[500 50 800 400]);
% sp(1) = subplot(121);
% plot(time2,mud(:,2),'-','linewidth',1.5,...
%     'Color','k'), hold on
% plot(time2,mud(:,1),'-','linewidth',1.5,...
%     'Color','r')
% sp(2) = subplot(122);
% pp(1) = plot(time2,fr(:,2),'-','linewidth',1.5,...
%     'Color','k'); hold on
% pp(2) = plot(time2,fr(:,1),'-','linewidth',1.5,...
%     'Color','r');
% leg = legend(pp,{'Lower';'Upper'});
% %plot adjustments
% set(sp(1),'position',[0.09 0.13 0.42 0.8])
% set(sp(2),'position',[0.55 0.13 0.42 0.8],...
%     'yticklabel',[])
% set(sp,'Xlim',[0 10],'ylim',[-0.5 0.4],...
%     'ytick',-0.5:0.25:0.5,'linewidth',1.5)
% set(leg,'position',[0.89 0.82 0.05 0.05])
% ylabel(sp(1),'Velocity (ms^-^1)')
% xlabel(sp(1),'Time (s)'),xlabel(sp(2),'Time (s)')
% title(sp(1),'Mudflat'),title(sp(2),'Mangrove Fringe')
% export_fig('d:\Projects\Documents\Writing\TOSpaper\Figures\Vel_Eps_SSC\Part1','-pdf','-nocrop')

f2 = figure(2);
cc = brewermap([],'*RdYlBu');
colormap(cc);
set(f2,'PaperOrientation','portrait',...
    'position',[500 50   800   1000]);
sp(1) = subplot(421);
imagesc(Stime(b(1):e(1)),rb,eps(:,b(1):e(1))),
set(gca,'ydir','normal','xlim',[Stime(b(1)) Stime(e(1))],...
    'xtick',Stime(b(1)):datenum(0,0,0,0,90,0):Stime(e(1)))
caxis([-5 -3.3])
datetick('x','HH:MM','keepticks','keeplimits')

sp(2) = subplot(422);
plot(time2(r(1):q(1)),SSC2(r(1):q(1)),...
    'linewidth',1.5,'color','k')
set(gca,'xlim',[time2(r(1)) time2(q(1))],...
    'xtick',time2(r(1)):datenum(0,0,0,0,90,0):time2(q(1)))
datetick('x','HH:MM','keepticks','keeplimits')

sp(3) = subplot(423);
imagesc(Stime(b(2):e(2)),rb,eps(:,b(2):e(2))),
set(gca,'ydir','normal','xlim',[Stime(b(2)) Stime(e(2))],...
    'xtick',Stime(b(2)):datenum(0,0,0,0,60,0):Stime(e(2)))
caxis([-5 -3.3])
datetick('x','HH:MM','keepticks','keeplimits')

sp(4) = subplot(424);
plot(time2(r(2):q(2)),SSC2(r(2):q(2)),...
    'linewidth',1.5,'color','k')
set(gca,'xlim',[time2(r(2)) time2(q(2))],...
    'xtick',time2(r(2)):datenum(0,0,0,0,60,0):time2(q(2)))
datetick('x','HH:MM','keepticks','keeplimits')

sp(5) = subplot(425);
imagesc(Stime(b(3):e(3)),rb,eps(:,b(3):e(3))),
set(gca,'ydir','normal','xlim',[Stime(b(3)) Stime(e(3))],...
    'xtick',Stime(b(3)):datenum(0,0,0,0,5,0):Stime(e(3)))
caxis([-5 -3.3])
datetick('x','HH:MM','keepticks','keeplimits')

sp(6) = subplot(426);
plot(time2(r(3):q(3)),SSC2(r(3):q(3)),...
    'linewidth',1.5,'color','k')
set(gca,'xlim',[time2(r(3)) time2(q(3))],...
    'xtick',time2(r(3)):datenum(0,0,0,0,5,0):time2(q(3)))
datetick('x','HH:MM','keepticks','keeplimits')

sp(7) = subplot(427);
imagesc(Stime(b(4):e(4)),rb,eps(:,b(4):e(4))),
set(gca,'ydir','normal','xlim',[Stime(b(4)) Stime(e(4))],...
    'xtick',Stime(b(4)):datenum(0,0,0,0,1,0):Stime(e(4)))
caxis([-5 -3.3])
datetick('x','HH:MM','keepticks','keeplimits')
cb = colorbar('location','southoutside');

sp(8) = subplot(428);
plot(time2(r(4):q(4)),SSC2(r(4):q(4)),...
    'linewidth',1.5,'color','k')
set(gca,'xlim',[time2(r(4)) time2(q(4))],...
    'xtick',time2(r(4)):datenum(0,0,0,0,1,0):time2(q(4)))
datetick('x','HH:MM','keepticks','keeplimits')


%Plot adjustments
set(sp(1),'position',[0.1 0.81 0.27 0.18])
set(sp(2),'position',[0.5 0.81 0.27 0.18],...
    'ytick',250:750:2000,'ylim',[0 2000])
set(sp(3),'position',[0.16 0.58 0.27 0.18])
set(sp(4),'position',[0.56 0.58 0.27 0.18],...
    'ytick',250:500:1250,'ylim',[0 1500])
set(sp(5),'position',[0.23 0.36 0.27 0.18])
set(sp(6),'position',[0.64 0.36 0.27 0.18],...
    'ytick',900:100:1100,'ylim',[850 1150])
set(sp(7),'position',[0.3 0.14 0.27 0.18])
set(sp(8),'position',[0.7 0.14 0.27 0.18],...
    'ytick',900:100:1000,'ylim',[850 1100])
set(cb,'position',[0.3 0.06 0.27 0.02],...
    'xtick',-5:.5:-3.5,'xlim',[-5 -3.5])
xlabel(cb,'log_1_0(\epsilon) (Wkg^-^1)')
ylabel(sp(1),'HAB (m)')
ylabel(sp(2),'SSC (mg l^-^1)')
ylabel(sp(3),'HAB (m)')
ylabel(sp(4),'SSC (mg l^-^1)')
ylabel(sp(5),'HAB (m)')
ylabel(sp(6),'SSC (mg l^-^1)')
ylabel(sp(7),'HAB (m)')
xlabel(sp(7),['Time on ' datestr(Stime(1),'dd-mm-yyyy')])
ylabel(sp(8),'SSC (mg l^-^1)')
xlabel(sp(8),['Time on ' datestr(Stime(1),'dd-mm-yyyy')])
prettyfigures('text',13,'labels',14,'box',1)
export_fig('d:\Projects\Documents\Writing\TOSpaper\Figures\Vel_Eps_SSC\Part2_v2','-pdf','-nocrop')





