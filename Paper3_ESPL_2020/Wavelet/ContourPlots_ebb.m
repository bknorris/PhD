clear,close all
load('G:\Mekong_W2015\DataAnalysis\Paper3\Wavelet\CombEventData_flood');
sfdir = 'f:\GradSchool\DataAnalysis\Paper3\WorkingFigures\Wavelet\';
units = {'min';'m';'m';'m/s';'m/s'};
cl = [207 176 126;
    60 166 74;
    4 76 41]./255;
symb = {'o';'^';'s'};
fn = fieldnames(dat);
%Following Joss's recommendation to make contours for color plots (as my
%former ones didn't work out so well!)
%% UMAG and UORB colored by BDMED
sp = zeros(3,1);
f1 = figure(1);
set(f1,'PaperOrientation','portrait',...
    'position',[400 100   1000   400],...
    'renderer','painters');    
factor = 35; %joss's npred
range = 1; %joss's perc_range
cmap = brewermap(factor,'RdBu');cmap = flipud(cmap);
for i = 1:3
    sp(i) = subplot(1,3,i);
    %create a patch for background color
    px = [0 0.5 0.5 0];
    py = [0.06 0.06 0.45 0.45];
    patch(px,py,'k'),hold on
    umag = [dat.(fn{i}).wave.umag; dat.(fn{i}).ig.umag];
    bdmed = [dat.(fn{i}).wave.bdmed; dat.(fn{i}).ig.bdmed];
    uorb = [dat.(fn{i}).wave.orbwv; dat.(fn{i}).ig.orbwv];
    %prep data for contour
    mx = linspace(min(umag)*range,max(umag)*range,factor);
    my = linspace(min(uorb)*range,max(uorb)*range,factor);
    [vx,vy] = meshgrid(mx,my);
    v = griddata(umag,uorb,bdmed,vx,vy,'cubic');
    p = contourf(vx,vy,v);
    colormap(cmap);
    % shading(gca,'interp')
    caxis([-0.03 0.03])
end
set(sp(1),'position',[0.1 0.18 0.2 0.72],...
    'xlim',[0.06 0.45],'ylim',[0.15 0.4])
set(sp(2),'position',[0.37 0.18 0.2 0.72],...
    'xlim',[0.06 0.45],'ylim',[0.15 0.4])
set(sp(3),'position',[0.64 0.18 0.2 0.72],...
    'xlim',[0.05 0.2],'ylim',[0.15 0.25],...
    'ytick',0.15:0.025:0.25)
title(sp(1),'Mudflat'),title(sp(2),'Fringe'),title(sp(3),'Forest')
xtext = sprintf(['Cross-shore velocity' '\n' 'magnitude [m/s]']);
xlabel(sp(1),xtext),xlabel(sp(2),xtext),xlabel(sp(3),xtext)
ylabel(sp(1),'Wave Orbital Velocity [m/s]')
cb = colorbar;
ytext = sprintf(['Median elev. relative' '\n' 'to event beginning [m]']);
ylabel(cb,ytext)
set(cb,'position',[0.86 0.18 0.02 0.72])
prettyfigures('text',12,'labels',13,'box',1)
% export_fig([sfdir 'F_contour_umagUorb'],'-png')
%%%
