%plot script for TKE output- SimpleStructureFunct.m
savedatdir = 'c:\Users\bkn5\Projects\Mekong_W2015\DataAnalysis\Turbulence\Paper1\';
figdir = 'c:\Users\bkn5\Projects\Mekong_W2015\Figures\Turbulence\Statistics\';

load([savedatdir 'F2F2day2TKE.mat'])

%%%% VP1
[m,n] = size(Stat.vpro1.beam1.A);
time = 10:10:10*n;

%TKE color plot
f1 = figure(1);
set(f1,'PaperOrientation','portrait',...
    'position',[200 200   1000   600]);
set(gcf,'color','w','PaperPositionMode','auto')
epsilon = (Stat.vpro1.beam1.E+Stat.vpro1.beam3.E)./2;
epserr = (Stat.vpro1.beam1.Emaxerr+Stat.vpro1.beam3.Emaxerr)./2;

p(1) = subplot(131);
imagesc(time,0:30,epsilon-epserr);
set(gca,'Xlim',[0 10*n],'YLim',[0 30])
title('Minimum \epsilon')
xlabel('Minutes Elapsed')
ylabel('Depth Bins (0.001m)')
caxis([0 5E-3])

p(2) = subplot(132);
imagesc(time,0:30,epsilon)
set(gca,'Xlim',[0 10*n],'YLim',[0 30],'YTickLabel',[])
title('\epsilon')
xlabel('Minutes Elapsed')
caxis([0 5E-3])

p(3) = subplot(133);
imagesc(time,0:30,epsilon+epserr)
set(gca,'Xlim',[0 10*n],'YLim',[0 30],'YTickLabel',[])
title('Maximum \epsilon')
xlabel('Minutes Elapsed')
set(p(1),'position',[0.1 0.22 0.25 0.7])
set(p(2),'position',[0.38 0.22 0.25 0.7])
set(p(3),'position',[0.66 0.22 0.25 0.7])
caxis([0 5E-3])
cb = colorbar('south');
set(cb,'position',[0.25 0.02, 0.54 0.05]),xlabel(cb,'W/kg')
suptitle('VP1 F2F2 Day 2 \epsilon Confidence Estimates')
export_fig([figdir 'F2F2day2VP1epsCI'],'-jpeg')

%Best fit predicted slope
bfit = (Stat.vpro1.beam1.bfit+Stat.vpro1.beam3.bfit)./2;
twothirds = (2/3)*ones(n,1);

f2 = figure(2);
set(f2,'PaperOrientation','portrait',...
    'position',[200 200   800  600]);
set(gcf,'color','w','PaperPositionMode','auto')
c = linspecer(n,'Sequential');
plot(time,twothirds,'--k','LineWidth',1.5), hold on
for i = 1:n
    plot(time,bfit(i,:),'o','Color',c(i,:),'MarkerFaceColor',c(i,:),'MarkerEdgeColor','k')
end
set(gca,'Xlim',[10 10*n],'YLim',[0.5 1.5])
xlabel('Minutes Elapsed')
ylabel('Slope')
title('VP1 F2F2 Day 2 Best-Fit Slope Estimates for r^2^/^3')
export_fig([figdir 'F2F2day2VP1bfslopes'],'-jpeg')

%r^2/3 r^2 values
RSQ = (Stat.vpro1.beam1.RSQ+Stat.vpro1.beam3.RSQ)./2;

f3 = figure(3);
set(f3,'PaperOrientation','portrait',...
    'position',[200 200   800  600]);
set(gcf,'color','w','PaperPositionMode','auto')
c = linspecer(n,'Sequential');
for i = 1:n
    plot(time,RSQ(i,:),'o','Color',c(i,:),'MarkerFaceColor',c(i,:),'MarkerEdgeColor','k'), hold on
end
set(gca,'Xlim',[10 10*n],'YLim',[0.95 1])
xlabel('Minutes Elapsed')
ylabel('r^2 value')
title('VP1 F2F2 Day 2 r^2 fits for r^2^/^3')
export_fig([figdir 'F2F2day2VP1rsq'],'-jpeg')

%example r^(2/3) fits versus r^b fits; plot only 5 bins
f4 = figure(4);
set(f4,'PaperOrientation','portrait',...
    'position',[200 200   800  600]);
set(gcf,'color','w','PaperPositionMode','auto')
c = autumn(5);
D = (Stat.vpro1.beam1.D(1:5,1:5)+Stat.vpro1.beam3.D(1:5,1:5))./2;
Yfit = (Stat.vpro1.beam1.Yfit(1:5,1:5)+Stat.vpro1.beam3.Yfit(1:5,1:5))./2;
% yfit = (Stat.vpro1.beam1.yfit(1:5,1:5)+Stat.vpro1.beam3.yfit(1:5,1:5))./2;
R = (Stat.vpro1.beam1.R(:,1)+Stat.vpro1.beam3.R(:,1))./2;
% r = (Stat.vpro1.beam1.r(:,1)+Stat.vpro1.beam3.r(:,1))./2;
for ii = 1:5
    plot(R,D(ii,:),'o','MarkerFaceColor',c(ii,:),'MarkerEdgeColor','k'), hold on
    plot(R,Yfit(ii,:),'Color',c(ii,:),'LineWidth',1.5)
%     plot(r,D(ii,:),'d','Color',c(ii,:)), hold on
%     plot(r,yfit(ii,:),'Color',c(ii,:))
end
xlabel('Along Beam Distance r^2^/^3')
ylabel('D(z,r)')
title('VP1 F2F2 Day 2 D(z,r) = N+Ar^2^/^3 fits')
export_fig([figdir 'F2F2day2VP1linfit'],'-jpeg')

