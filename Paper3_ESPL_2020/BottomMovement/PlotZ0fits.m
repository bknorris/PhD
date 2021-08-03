clear, close all
load('d:\Mekong_W2015\DataAnalysis\Paper3\BedStress\06-03-15\F2F2_06_vpro1.mat')
figure
plot(dat.z0)
idx(1) = 210; %a good fit of umag, log(zuv)
idx(2) = 3999; %a poor fit of umag, log(zuv)
f1 = figure(1);
set(f1,'PaperOrientation','portrait',...
    'position',[400 100   900   400],...
    'renderer','painters');hold on
sp = zeros(1,2);
for i = 1:2
    sp(i) = subplot(1,2,i);
    zuv = dat.z;
    zid = dat.z<=0.005;
    zuv(zid) = [];
    umag = dat.umag(idx(i),:);
    umag(zid) = [];
    xs = log(zuv);ys = umag;
    pf = polyfit(xs,ys,1);
    pv = polyval(pf,xs);
    plot(xs,ys,'ok','linewidth',1.5),hold on
    plot(xs,pv,'r-','linewidth',1.5)
    eqn = sprintf('u(z) = %0.1dln(z)-%0.1d',pf(1),pf(2));
    tt(1) = text(0.1,0.95,eqn,'units','normalized');
    z0 = exp(-pf(2)/pf(1));
    ztxt = sprintf('Estimated z0: %0.2d m',z0);
    tt(2) = text(0.1,0.88,ztxt,'units','normalized');
    set(tt,'fontsize',11,'fontweight','bold')
end
set(sp(1),'position',[0.12 0.15 0.35 0.75])
set(sp(2),'position',[0.55 0.15 0.35 0.75])
title(sp(1),'A Good Fit')
title(sp(2),'A Poor Fit')
xlabel(sp(1),'ln(z) [m]')
xlabel(sp(2),'ln(z) [m]')
ylabel(sp(1),'u(z) [m/s]')
% prettyfigures('text',12,'labels',13,'box',1)
% sfdir = 'd:\Projects\Mekong_W2015\Figures\Paper3\BottomTracking\';
% export_fig([sfdir 'z0_fits_good_bad'],'-png')
