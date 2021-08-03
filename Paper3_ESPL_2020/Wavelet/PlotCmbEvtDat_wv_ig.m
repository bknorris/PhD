%Plot the combined event data
clear,close all
load('f:\Mekong_W2015\DataAnalysis\Paper3\Wavelet\CombEventData_flood.mat');
stage = 'flood';
sfdir = 'g:\GradSchool\DataAnalysis\Paper3\WorkingFigures\Wavelet\';
units = {'min';'m';'m';'m/s';'m/s'};
cl = [207 176 126;
    60 166 74;
    4 76 41]./255;
symb = {'o';'^';'s'};
fn = fieldnames(dat);
%% Median Bed Level
%Event Level
figure
eb = zeros(3,1);
for i = 1:length(fn)
    eventl = [dat.(fn{i}).wave.eventl; dat.(fn{i}).ig.eventl];
    bins = linspace(min(eventl),max(eventl),10);
    bdmed = [dat.(fn{i}).wave.bdmed; dat.(fn{i}).ig.bdmed];
    [b,~,q1,q3] = binmedian(eventl,bdmed,bins);
    eb(i) = errorbar(bins,b,q1,q3,symb{i},...
        'color',cl(i,:),...
        'linewidth',1.5,...
        'markersize',8,...
        'markeredgecolor','k',...
        'markerfacecolor',cl(i,:));hold on
    
end
xlabel('Event Length [min]')
ytext = sprintf(['Median elev. relative' '\n' 'to event beginning [m]']);
ylabel(ytext)
leg = legend(eb,fn);
set(leg,'position',[0.27 0.2 0.05 0.05])
prettyfigures('text',12,'labels',13,'box',1)
export_fig([sfdir 'CmbData_' stage '_eventl'],'-png')
%%%
%Wave Orbital Velocity
figure
eb = zeros(3,1);
for i = 1:length(fn)
    uorb = [dat.(fn{i}).wave.orbwv; dat.(fn{i}).ig.orbwv];
    bins = linspace(min(uorb),max(uorb),10);
    bdmed = [dat.(fn{i}).wave.bdmed; dat.(fn{i}).ig.bdmed];
    [b,~,q1,q3] = binmedian(uorb,bdmed,bins);
    eb(i) = errorbar(bins,b,q1,q3,symb{i},...
        'color',cl(i,:),...
        'linewidth',1.5,...
        'markersize',8,...
        'markeredgecolor','k',...
        'markerfacecolor',cl(i,:));hold on
end
xlabel('Orbital Wave Velocity [m/s]')
ytext = sprintf(['Median elev. relative' '\n' 'to event beginning [m]']);
ylabel(ytext)
leg = legend(eb,fn);
set(leg,'position',[0.25 0.2 0.05 0.05])
prettyfigures('text',12,'labels',13,'box',1)
export_fig([sfdir 'CmbData_' stage '_uorb'],'-png')
%%%
%Cross-shore velocity magnitude
figure
eb = zeros(3,1);
for i = 1:length(fn)
    umag = [dat.(fn{i}).wave.umag; dat.(fn{i}).ig.umag];
    bins = linspace(min(umag),max(umag),10);
    bdmed = [dat.(fn{i}).wave.bdmed; dat.(fn{i}).ig.bdmed];
    [b,~,q1,q3] = binmedian(umag,bdmed,bins);
    eb(i) = errorbar(bins,b,q1,q3,symb{i},...
        'color',cl(i,:),...
        'linewidth',1.5,...
        'markersize',8,...
        'markeredgecolor','k',...
        'markerfacecolor',cl(i,:));hold on
end
xlabel('Cross-shore velocity magnitude [m/s]')
ytext = sprintf(['Median elev. relative' '\n' 'to event beginning [m]']);
ylabel(ytext)
leg = legend(eb,fn);
set(leg,'position',[0.25 0.2 0.05 0.05])
prettyfigures('text',12,'labels',13,'box',1)
export_fig([sfdir 'CmbData_' stage '_umag'],'-png')
%%%
%Significant Wave Height
figure
eb = zeros(3,1);
for i = 1:length(fn)
    Hs = [dat.(fn{i}).wave.sigh; dat.(fn{i}).ig.sigh];
    bins = linspace(min(Hs),max(Hs),10);
    bdmed = [dat.(fn{i}).wave.bdmed; dat.(fn{i}).ig.bdmed];
    [b,~,q1,q3] = binmedian(Hs,bdmed,bins);
    eb(i) = errorbar(bins,b,q1,q3,symb{i},...
        'color',cl(i,:),...
        'linewidth',1.5,...
        'markersize',8,...
        'markeredgecolor','k',...
        'markerfacecolor',cl(i,:));hold on
end
xlabel('H_s [m]')
ytext = sprintf(['Median elev. relative' '\n' 'to event beginning [m]']);
ylabel(ytext)
leg = legend(eb,fn);
set(leg,'position',[0.25 0.2 0.05 0.05])
prettyfigures('text',12,'labels',13,'box',1)
export_fig([sfdir 'CmbData_' stage '_Hs'],'-png')
%%%
%Depth
figure
eb = zeros(3,1);
for i = 1:length(fn)
    H = [dat.(fn{i}).wave.depth; dat.(fn{i}).ig.depth];
    bins = linspace(min(H),max(H),10);
    bdmed = [dat.(fn{i}).wave.bdmed; dat.(fn{i}).ig.bdmed];
    [b,~,q1,q3] = binmedian(H,bdmed,bins);
    eb(i) = errorbar(bins,b,q1,q3,symb{i},...
        'color',cl(i,:),...
        'linewidth',1.5,...
        'markersize',8,...
        'markeredgecolor','k',...
        'markerfacecolor',cl(i,:));hold on
end
xlabel('Water Depth [m]')
ytext = sprintf(['Median elev. relative' '\n' 'to event beginning [m]']);
ylabel(ytext)
leg = legend(eb,fn);
set(leg,'position',[0.25 0.2 0.05 0.05])
prettyfigures('text',12,'labels',13,'box',1)
export_fig([sfdir 'CmbData_' stage '_depth'],'-png')
%%%
%Tau_b
figure
eb = zeros(3,1);
for i = 1:length(fn)
    taub = [dat.(fn{i}).wave.taub; dat.(fn{i}).ig.taub];
    bins = linspace(min(taub),max(taub),10);
    bdmed = [dat.(fn{i}).wave.bdmed; dat.(fn{i}).ig.bdmed];
    [b,~,q1,q3] = binmedian(taub,bdmed,bins);
    eb(i) = errorbar(bins,b,q1,q3,symb{i},...
        'color',cl(i,:),...
        'linewidth',1.5,...
        'markersize',8,...
        'markeredgecolor','k',...
        'markerfacecolor',cl(i,:));hold on
end
xlabel('Bed Shear Stress [Pa]')
ytext = sprintf(['Median elev. relative' '\n' 'to event beginning [m]']);
ylabel(ytext)
leg = legend(eb,fn);
set(leg,'position',[0.25 0.2 0.05 0.05])
prettyfigures('text',12,'labels',13,'box',1)
export_fig([sfdir 'CmbData_' stage '_taub'],'-png')
%%%
%epsilon
figure
eb = zeros(3,1);
for i = 1:length(fn)
    eps = [dat.(fn{i}).wave.eps; dat.(fn{i}).ig.eps];
    bins = linspace(min(eps),max(eps),10);
    bdmed = [dat.(fn{i}).wave.bdmed; dat.(fn{i}).ig.bdmed];
    [b,~,q1,q3] = binmedian(eps,bdmed,bins);
    eb(i) = errorbar(bins,b,q1,q3,symb{i},...
        'color',cl(i,:),...
        'linewidth',1.5,...
        'markersize',8,...
        'markeredgecolor','k',...
        'markerfacecolor',cl(i,:));hold on
end
xlabel('TKE Dissipation Rate [m^2/s^2]')
ytext = sprintf(['Median elev. relative' '\n' 'to event beginning [m]']);
ylabel(ytext)
leg = legend(eb,fn);
set(leg,'position',[0.25 0.2 0.05 0.05])
prettyfigures('text',12,'labels',13,'box',1)
export_fig([sfdir 'CmbData_' stage '_eps'],'-png')
%%%
%usquared
figure
eb = zeros(3,1);
for i = 1:length(fn)
    usqd = [dat.(fn{i}).wave.usqd; dat.(fn{i}).ig.usqd];
    bins = linspace(min(usqd),max(usqd),10);
    bdmed = [dat.(fn{i}).wave.bdmed; dat.(fn{i}).ig.bdmed];
    [b,~,q1,q3] = binmedian(usqd,bdmed,bins);
    eb(i) = errorbar(bins,b,q1,q3,symb{i},...
        'color',cl(i,:),...
        'linewidth',1.5,...
        'markersize',8,...
        'markeredgecolor','k',...
        'markerfacecolor',cl(i,:));hold on
end
xlabel('Velocity Squared [m^2/s^2]')
ytext = sprintf(['Median elev. relative' '\n' 'to event beginning [m]']);
ylabel(ytext)
leg = legend(eb,fn);
set(leg,'position',[0.25 0.2 0.05 0.05])
prettyfigures('text',12,'labels',13,'box',1)
export_fig([sfdir 'CmbData_' stage '_usqd'],'-png')
%% Normalized Bed Variance
%Event Level
figure
eb = zeros(3,1);
for i = 1:length(fn)
    dfn = fieldnames(dat.(fn{i}));
    eventl = [dat.(fn{i}).wave.eventl; dat.(fn{i}).ig.eventl];
    bins = linspace(min(eventl),max(eventl),10);
    bdvar = [dat.(fn{i}).wave.bdvar; dat.(fn{i}).ig.bdvar];
    [b,~,q1,q3] = binmedian(eventl,bdvar,bins);
    eb(i) = errorbar(bins,b,q1,q3,symb{i},...
        'color',cl(i,:),...
        'linewidth',1.5,...
        'markersize',8,...
        'markeredgecolor','k',...
        'markerfacecolor',cl(i,:));hold on
    
end
set(gca,'yscale','log')
%%%
%Wave Orbital Velocity
figure
eb = zeros(3,1);
for i = 1:length(fn)
    uorb = [dat.(fn{i}).wave.orbwv; dat.(fn{i}).ig.orbwv];
    bins = linspace(min(uorb),max(uorb),10);
    bdvar = [dat.(fn{i}).wave.bdvar; dat.(fn{i}).ig.bdvar];
    [b,~,q1,q3] = binmedian(uorb,bdvar,bins);
    eb(i) = errorbar(bins,b,q1,q3,symb{i},...
        'color',cl(i,:),...
        'linewidth',1.5,...
        'markersize',8,...
        'markeredgecolor','k',...
        'markerfacecolor',cl(i,:));hold on
end
set(gca,'yscale','log')
%%%
%Cross-shore velocity magnitude
figure
eb = zeros(3,1);
for i = 1:length(fn)
    umag = [dat.(fn{i}).wave.umag; dat.(fn{i}).ig.umag];
    bins = linspace(min(uorb),max(uorb),10);
    bdvar = [dat.(fn{i}).wave.bdvar; dat.(fn{i}).ig.bdvar];
    [b,~,q1,q3] = binmedian(umag,bdvar,bins);
    eb(i) = errorbar(bins,b,q1,q3,symb{i},...
        'color',cl(i,:),...
        'linewidth',1.5,...
        'markersize',8,...
        'markeredgecolor','k',...
        'markerfacecolor',cl(i,:));hold on
end
set(gca,'yscale','log')
%% Phase vs. Event Length
figure
eb = zeros(3,1);
for i = 1:length(fn)
    eventl = [dat.(fn{i}).wave.eventl; dat.(fn{i}).ig.eventl];
    bins = linspace(min(eventl),max(eventl),10);
    phase = [dat.(fn{i}).wave.phase; dat.(fn{i}).ig.phase];
    phaser = deg2rad(phase);
    phase = exp(sqrt(-1).*phaser);phase = rad2deg(phase);
    [b,~,q1,q3] = binmedian(eventl,phase,bins);
    eb(i) = errorbar(bins,b,b-q1,b-q3,symb{i},...
        'color',cl(i,:),...
        'linewidth',1.5,...
        'markersize',8,...
        'markeredgecolor','k',...
        'markerfacecolor',cl(i,:));hold on
end
set(gca,'ylim',[-180 180],'ytick',-180:90:180)
xlabel('Event Length [min]')
ylabel('Phase of Wavelet [deg]')
leg = legend(eb,fn);
set(leg,'position',[0.78 0.8 0.05 0.05])
prettyfigures('text',12,'labels',13,'box',1)
export_fig([sfdir 'CmbData_' stage '_phase'],'-png')
%% Contour Plots
%Following Joss's recommendation to make contours for color plots (as my
%former ones didn't work out so well!)
%USQD and UORB colored by BDMED
sp1 = zeros(3,1);
sp2 = zeros(3,1);
sp3 = zeros(3,1);
f1 = figure(1);
set(f1,'PaperOrientation','portrait',...
    'position',[400 100   1000   1000],...
    'renderer','opengl');
factor = 20; %joss's npred    
range = 1; %joss's perc_range
cmap = brewermap(18,'RdBu');cmap = flipud(cmap);
pl = [1 2 3;4 5 6;7 8 9];
for i = 1:3
    sp1(i) = subplot(3,3,pl(1,i));
    %create a patch for background color
%     px = [0 1.5 1.5 0];
%     py = [0 0 0.7 0.7];
%     patch(px,py,'k'),hold on
    depth = [dat.(fn{i}).wave.depth; dat.(fn{i}).ig.depth];
    eps = [dat.(fn{i}).wave.eps; dat.(fn{i}).ig.eps];
    bdmed = [dat.(fn{i}).wave.bdmed; dat.(fn{i}).ig.bdmed];
    %prep data for contour
    mx = linspace(min(eps)*range,max(eps)*range,factor);
    my = linspace(min(depth)*range,max(depth)*range,factor);
    [vx,vy] = meshgrid(mx,my);
    v = griddata(eps,depth,bdmed,vx,vy,'cubic');
    contourf(vx,vy,v);
    colormap(cmap);
    caxis([-0.03 0.03])
    %%%
    sp2(i) = subplot(3,3,pl(2,i));
    %create a patch for background color
%     px = [0 1.5 1.5 0];
%     py = [0 0 0.7 0.7];
%     patch(px,py,'k'),hold on
    depth = [dat.(fn{i}).wave.depth; dat.(fn{i}).ig.depth];
    taub = [dat.(fn{i}).wave.taub; dat.(fn{i}).ig.taub];
    bdmed = [dat.(fn{i}).wave.bdmed; dat.(fn{i}).ig.bdmed];
    %prep data for contour
    mx = linspace(min(taub)*range,max(taub)*range,factor);
    my = linspace(min(depth)*range,max(depth)*range,factor);
    [vx,vy] = meshgrid(mx,my);
    v = griddata(taub,depth,bdmed,vx,vy,'cubic');
    contourf(vx,vy,v);
    colormap(cmap);
    caxis([-0.03 0.03])
    %%%
    sp3(i) = subplot(3,3,pl(3,i));
    %create a patch for background color
%     px = [0 1.5 1.5 0];
%     py = [0 0 0.7 0.7];
%     patch(px,py,'k'),hold on
    usqd = [dat.(fn{i}).wave.usqd; dat.(fn{i}).ig.usqd];
    uorb = [dat.(fn{i}).wave.orbwv; dat.(fn{i}).ig.orbwv];
    bdmed = [dat.(fn{i}).wave.bdmed; dat.(fn{i}).ig.bdmed];
    %prep data for contour
    mx = linspace(min(usqd)*range,max(usqd)*range,factor);
    my = linspace(min(uorb)*range,max(uorb)*range,factor);
    [vx,vy] = meshgrid(mx,my);
    v = griddata(usqd,uorb,bdmed,vx,vy,'cubic');
    contourf(vx,vy,v);
    colormap(cmap);
    caxis([-0.03 0.03])
end
%Top Row - Epsilon & Depth
set(sp1(1),'position',[0.1 0.74 0.24 0.22],...
    'xlim',[0 6E-4],'ylim',[0.2 1.2],...
    'ytick',0:0.2:1.2)
set(sp1(2),'position',[0.4 0.74 0.24 0.22],...
    'xlim',[0 6E-4],'ylim',[0.2 1.2],...
    'ytick',0:0.2:1.2)
set(sp1(3),'position',[0.7 0.74 0.24 0.22],...
    'xlim',[0 1E-4],'ylim',[0.2 0.8],...
    'ytick',0:0.2:1.2)
xlabel(sp1(2),'\epsilon [Wkg^{-1}]')
ylabel(sp1(1),'Water Depth [m]')
%Middle Row - Taub & Depth
set(sp2(1),'position',[0.1 0.45 0.24 0.22],...
    'xlim',[0.3 1.3],'ylim',[0.2 1.2],...
    'ytick',0:0.2:1.2)
set(sp2(2),'position',[0.4 0.45 0.24 0.22],...
    'xlim',[0.3 1.3],'ylim',[0.2 1.2],...
    'ytick',0:0.2:1.2)
set(sp2(3),'position',[0.7 0.45 0.24 0.22],...
    'xlim',[0 0.45],'ylim',[0.2 0.8],...
    'ytick',0:0.2:1.2)
xlabel(sp2(2),'\tau_{b} [Pa]')
ylabel(sp2(1),'Water Depth [m]')
%Bottom Row - Usqd & uorb
set(sp3(1),'position',[0.1 0.15 0.24 0.23],...
    'xlim',[0.02 0.25],'ylim',[0.05 0.35],...
    'ytick',0.05:0.1:0.35)
set(sp3(2),'position',[0.4 0.15 0.24 0.23],...
    'xlim',[0.02 0.25],'ylim',[0.05 0.35],...
    'ytick',0.05:0.1:0.35)
set(sp3(3),'position',[0.7 0.15 0.24 0.23],...
    'xlim',[0.005 0.045],'ylim',[0 0.25],...
    'ytick',0:0.1:25)
xlabel(sp3(2),'\tau_{b} [Pa]')
ylabel(sp3(1),'u^2 [m^2s^{-2}]')
%%%
title(sp1(1),'Mudflat'),title(sp1(2),'Fringe'),title(sp1(3),'Forest')
cb = colorbar('southoutside');
set(cb,'position',[0.1 0.06 0.84 0.02])
xtext = 'Median bed elevation relative to event beginning [m]';
xlabel(cb,xtext)
prettyfigures('text',12,'labels',13,'box',1)
export_fig([sfdir 'Contour_' stage '_HepsTaubUsqd'],'-pdf')

