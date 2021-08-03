clear,close all
load('d:\Mekong_W2015\DataAnalysis\Paper3\Wavelet\CombEventData_flood.mat');
sfdir = 'e:\GradSchool\DataAnalysis\Paper3\WorkingFigures\Wavelet\';
units = {'min';'m';'m';'m/s';'m/s'};
cl = [207 176 126;
    60 166 74;
    4 76 41]./255;
symb = {'o';'^';'s'};
fn = fieldnames(dat);
%Following Joss's recommendation to make contours for color plots (as my
%former ones didn't work out so well!)
%% Eps, H and deltbd
sp = zeros(3,2);
eb = zeros(3,2);
f2 = figure(2);
set(f2,'PaperOrientation','portrait',...
    'position',[400 100   1000   500],...
    'renderer','painters');
w = [1 4; 2 5; 3 6];
factor = 30; %joss's npred
range = 1; %joss's perc_range
cmap = brewermap(factor,'RdBu');cmap = flipud(cmap);
for i = 1:3
    %Accretion
    sp(w(i,1)) = subplot(2,3,w(i,1));
    %create a patch for background color
%     px = [0 0.5 0.5 0];
%     py = [0.06 0.06 0.45 0.45];
%     patch(px,py,'k'),hold on
    eps = [dat.(fn{i}).wave.eps; dat.(fn{i}).ig.eps];
    H = [dat.(fn{i}).wave.depth; dat.(fn{i}).ig.depth];
    deltbd = [dat.(fn{i}).wave.deltbd; dat.(fn{i}).ig.deltbd];
    id = find(deltbd>0.0015);
    x = eps(id);y = H(id);z = deltbd(id);
    %prep data for contour
    mx = linspace(min(x)*range,max(x)*range,factor);
    my = linspace(min(y)*range,max(y)*range,factor);
    [vx,vy] = meshgrid(mx,my);
    v = griddata(x,y,z,vx,vy,'cubic');
    contourf(vx,vy,v);
    colormap(cmap);
    % shading(gca,'interp')
    caxis([-0.03 0.03])
    
    %Erosion
    sp(w(i,2)) = subplot(2,3,w(i,2));
    %create a patch for background color
%     px = [0 0.5 0.5 0];
%     py = [0.06 0.06 0.45 0.45];
%     patch(px,py,'k'),hold on
    id = find(deltbd<-0.0015);
    x = eps(id);y = H(id);z = deltbd(id);
    %prep data for contour
    mx = linspace(min(x)*range,max(x)*range,factor);
    my = linspace(min(y)*range,max(y)*range,factor);
    [vx,vy] = meshgrid(mx,my);
    v = griddata(x,y,z,vx,vy,'cubic');
    contourf(vx,vy,v);
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
