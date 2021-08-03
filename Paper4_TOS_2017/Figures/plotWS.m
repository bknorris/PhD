%first, get bubble control

clear
close all
datdir = 'c:\Users\bkn5\Projects\Mekong_W2015\DataAnalysis\Spectra\Paper1\';
hta = load([datdir 'WaveStatsV5108_HTA.mat']);
vta = load([datdir 'WaveStatsV5109_VTA.mat']);
datdir = 'c:\Users\bkn5\Projects\Mekong_W2015\DataAnalysis\Paper1\';
figdir = 'c:\Users\bkn5\Projects\Mekong_W2015\Figures\Environment\';
times = load([datdir 'HTAtimes.mat']);
hta.times = times.HTA.times; clear HTA times
times = load([datdir 'VTAtimes.mat']);
vta.times = times.VTA.times; clear VTA times
e1 = hta.wvstats.time(end);t2 = vta.wvstats.time(1);int = datenum(0,0,0,0,2,0); %2 minute sampling gap
tgap = (e1:int:t2)';
dgap = NaN(length(tgap),1);
data.times = [hta.wvstats.time; tgap; vta.wvstats.time];
data.depth = [hta.wvstats.depth-0.18; dgap; vta.wvstats.depth-0.2];
data.hrmsp = [hta.wvstats.hrmsp; dgap; vta.wvstats.hrmsp];
data.udir = [hta.wvstats.udir; dgap; vta.wvstats.udir];
data.azm = [hta.wvstats.azm; dgap; vta.wvstats.azm];
data.u = [hta.wvstats.u; dgap; vta.wvstats.u];
data.v = [hta.wvstats.v; dgap; vta.wvstats.v];
%extract times of beginning-end of experiments to plot transect lines on
%the azimuth subplot, and vertical lines on all subplots denoting the
%transition between study sites
ex1b = data.times(1);ex1e = hta.wvstats.time(end);ex2b = vta.wvstats.time(1);ex2e = data.times(end);
vb.xs1 = repmat(ex1e,1,400);vb.ys = 0:399;vb.xs2 = repmat(ex2b,1,400);
hb.xs1 = linspace(ex1b,ex1e,length(hta.wvstats.time));hb.ys1 = repmat(340,1,length(hta.wvstats.time));
hb.xs2 = linspace(ex2b,ex2e,length(vta.wvstats.time));hb.ys2 = repmat(290,1,length(vta.wvstats.time));


f1 = figure(1);
set(f1,'PaperOrientation','portrait',...
    'position',[400 200   1000   800]);
set(gcf,'color','w','PaperPositionMode','auto')

s(1) = subplot(311);hold on
p(1) = patch([hta.times.t1 hta.times.e1 hta.times.e1 hta.times.t1],[0 0 2 2],[.8 .8 .8]);
p(2) = patch([hta.times.t2 hta.times.e2 hta.times.e2 hta.times.t2],[0 0 2 2],[.8 .8 .8]);
p(3) = patch([hta.times.t3 hta.times.e3 hta.times.e3 hta.times.t3],[0 0 2 2],[.8 .8 .8]);
p(4) = patch([vta.times.t2 vta.times.e2 vta.times.e2 vta.times.t2],[0 0 2 2],[.8 .8 .8]);
l(1) = line(vb.xs1,vb.ys,'LineWidth',2,'Color','k');
l(2) = line(vb.xs2,vb.ys,'LineWidth',2,'Color','k');
g(1) = plot(data.times,data.depth,'-','Color','k','LineWidth',1.5);hold off
ylabel('\bf\itDepth above Bed (m)','FontSize',14)
set(gca,'XTickLabel',[],'LineWidth',1.5,'FontSize',14,'XLim',[data.times(1) data.times(end)],...
    'YLim',[0 2],'YTick',0:0.5:2,'box','on')

s(2) = subplot(312);hold on
p(5) = patch([hta.times.t1 hta.times.e1 hta.times.e1 hta.times.t1],[0 0 2 2],[.8 .8 .8]);
p(6) = patch([hta.times.t2 hta.times.e2 hta.times.e2 hta.times.t2],[0 0 2 2],[.8 .8 .8]);
p(7) = patch([hta.times.t3 hta.times.e3 hta.times.e3 hta.times.t3],[0 0 2 2],[.8 .8 .8]);
p(8) = patch([vta.times.t2 vta.times.e2 vta.times.e2 vta.times.t2],[0 0 2 2],[.8 .8 .8]);
l(1) = line(vb.xs1,vb.ys,'LineWidth',1.5,'Color','k');
l(2) = line(vb.xs2,vb.ys,'LineWidth',1.5,'Color','k');
g(2) = plot(data.times,data.hrmsp,'-','Color','k','LineWidth',1.5,'Marker','.','MarkerSize',3);hold off
ylabel('\bf\itHrms (m)','FontSize',14)
set(gca,'XTickLabel',[],'LineWidth',1.5,'FontSize',14,'XLim',[data.times(1) data.times(end)],...
    'YLim',[0 2],'YTick',0:0.5:1.5,'box','on')


s(3) = subplot(313);hold on
p(9) = patch([hta.times.t1 hta.times.e1 hta.times.e1 hta.times.t1],[0 0 400 400],[.8 .8 .8]);
p(10) = patch([hta.times.t2 hta.times.e2 hta.times.e2 hta.times.t2],[0 0 400 400],[.8 .8 .8]);
p(11) = patch([hta.times.t3 hta.times.e3 hta.times.e3 hta.times.t3],[0 0 400 400],[.8 .8 .8]);
p(12) = patch([vta.times.t2 vta.times.e2 vta.times.e2 vta.times.t2],[0 0 400 400],[.8 .8 .8]);
l(1) = line(vb.xs1,vb.ys,'LineWidth',1.5,'Color','k');
l(2) = line(vb.xs2,vb.ys,'LineWidth',1.5,'Color','k');
g(3) = plot(data.times,data.udir,'o','Color','k','MarkerSize',4,'MarkerFaceColor','k');
g(4) = plot(data.times,(360-data.azm),'+','Color',[0.65 0.65 0.65],'MarkerSize',5,'MarkerFaceColor',[0.65 0.65 0.65]);
l(3) = line(hb.xs1,hb.ys1,'LineWidth',2,'Color','r');
l(4) = line(hb.xs2,hb.ys2,'LineWidth',2,'Color','r');hold off
set(p,'EdgeColor','none')
ylabel('\bf\itF_m (Wm^-^1)','FontSize',14)
set(gca,'LineWidth',1.5,'FontSize',14,'XLim',[data.times(1) data.times(end)],'YLim',[0 400],'YTick',0:90:360,'box','on')
xlabel('\bf\itDays in 03/2015','FontSize',14)
ylabel('\bf\itGeographic Direction (deg)','FontSize',14)
datetick('x','dd HH:MM','keepticks','keeplimits')


set(s(1),'position',[0.1 0.66 0.82 0.28])
set(s(2),'position',[0.1 0.38 0.82 0.28])
set(s(3),'position',[0.1 0.1 0.82 0.28])

fpath = figdir;fname = 'Wavestats_fullts';
% export_fig([fpath fname],'-png')
% disp(['Figure ' fname '.png saved'])

clear p l s 
%Figure 2
f2 = figure(2);
set(f2,'PaperOrientation','portrait',...
    'position',[400 200   1000   800]);
set(gcf,'color','w','PaperPositionMode','auto')

%%%%depth
s(1) = subplot(311);hold on
p(1) = patch([hta.times.t1 hta.times.e1 hta.times.e1 hta.times.t1],[0 0 2 2],[.8 .8 .8]);
p(2) = patch([hta.times.t2 hta.times.e2 hta.times.e2 hta.times.t2],[0 0 2 2],[.8 .8 .8]);
p(3) = patch([hta.times.t3 hta.times.e3 hta.times.e3 hta.times.t3],[0 0 2 2],[.8 .8 .8]);
p(4) = patch([vta.times.t2 vta.times.e2 vta.times.e2 vta.times.t2],[0 0 2 2],[.8 .8 .8]);
g(1) = plot(data.times,data.depth,'-','Color','k','LineWidth',2);hold off
%adjustments
ylabel('Depth (m)','FontSize',14,'FontName','Cambria')
datetick('x','dd','keepticks','keeplimits')
set(p,'EdgeColor','none')
set(gca,'LineWidth',1.5,'FontSize',14,'FontName','Cambria','XLim',[data.times(1) data.times(end)],...
    'YLim',[0 2],'YTick',0:0.5:2,'box','on')
%%%%wave height
s(2) = subplot(312);hold on
p(1) = patch([hta.times.t1 hta.times.e1 hta.times.e1 hta.times.t1],[0 0 2 2],[.8 .8 .8]);
p(2) = patch([hta.times.t2 hta.times.e2 hta.times.e2 hta.times.t2],[0 0 2 2],[.8 .8 .8]);
p(3) = patch([hta.times.t3 hta.times.e3 hta.times.e3 hta.times.t3],[0 0 2 2],[.8 .8 .8]);
p(4) = patch([vta.times.t2 vta.times.e2 vta.times.e2 vta.times.t2],[0 0 2 2],[.8 .8 .8]);
Hs = data.hrmsp*(sqrt(2));id = find(Hs < 0.01);Hs(id) = NaN;
g(1) = plot(data.times(1:10:length(data.times)),Hs(1:10:length(Hs)),'ok',...
    'MarkerSize',6,'MarkerFaceColor',[1 1 1],'MarkerEdgeColor',[0 0 0],...
    'LineWidth',1.5);hold off
%adjustments
ylabel('H_s (m)','FontSize',14,'FontName','Cambria')
datetick('x','dd','keepticks','keeplimits')
set(p,'EdgeColor','none')
set(gca,'LineWidth',1.5,'FontSize',14,'FontName','Cambria','XLim',[data.times(1) data.times(end)],...
    'YLim',[0 2],'YTick',0:0.5:2,'box','on')
%%%%wave dir
s(3) = subplot(313);hold on
p(1) = patch([hta.times.t1 hta.times.e1 hta.times.e1 hta.times.t1],[0 0 360 360],[.8 .8 .8]);
p(2) = patch([hta.times.t2 hta.times.e2 hta.times.e2 hta.times.t2],[0 0 360 360],[.8 .8 .8]);
p(3) = patch([hta.times.t3 hta.times.e3 hta.times.e3 hta.times.t3],[0 0 360 360],[.8 .8 .8]);
p(4) = patch([vta.times.t2 vta.times.e2 vta.times.e2 vta.times.t2],[0 0 360 360],[.8 .8 .8]);
l(1) = line(hb.xs1,hb.ys1,'LineStyle','--','LineWidth',2,'Color','k');
l(2) = line(hb.xs2,hb.ys2,'LineStyle','--','LineWidth',2,'Color','k');
g(1) = plot(data.times(1:20:length(data.times)),360-data.azm(1:20:length(data.azm)),'ok',...
    'MarkerSize',6,'MarkerFaceColor',[1 1 1],'MarkerEdgeColor',[0 0 0],...
    'LineWidth',1.5);hold off
%adjustments
ylabel('Wave Direction (Deg)','FontSize',14,'FontName','Cambria')
xlabel('March 2015','FontSize',14,'FontName','Cambria')
datetick('x','dd','keepticks','keeplimits')
set(p,'EdgeColor','none')
set(gca,'LineWidth',1.5,'FontSize',14,'FontName','Cambria','XLim',[data.times(1) data.times(end)],...
    'YLim',[180 360],'YTick',180:90:360,'box','on')

set(s(1),'position',[0.1 0.7 0.85 0.22])
set(s(2),'position',[0.1 0.4 0.85 0.22])
set(s(3),'position',[0.1 0.1 0.85 0.22])

fpath = figdir;fname = 'Wavestats_simple';
export_fig([fpath fname],'-jpeg','-nocrop')
disp(['Figure ' fname '.jpeg saved'])

% %%%%current direction
% s(3) = subplot(414);
% %plot current whisks
% xlim = [data.times(1) data.times(end)];
% ylim = [-0.3 0.3];
% nw = length(data.u);
% scale = 1;ps = 1;
% usc = scale.*data.u;
% vsc = scale.*data.v;
% pl = ps.*data.times;
% 
% hold on
% for id = 1:15:nw;
% plot([pl(id); usc(id)+pl(id)],[0; vsc(id)],...
%     'color','k','linewidth',1.5);
% end
% 
% min_x = xlim(1);
% max_x = xlim(2);
% min_y = ylim(1);
% max_y = ylim(2);
% 
% axis equal
% axis([min_x*ps max_x*ps min_y*scale max_y*scale])
% l = line([min_x*ps max_x*ps],[0 0],'Color','k','LineWidth',1.5);

