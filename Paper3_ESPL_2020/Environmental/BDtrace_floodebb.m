%Plot the bottom tracking for the experiment on 09-03-15, GM 2018 Norris et
%al.
clear,close all
ddir = 'e:\Mekong_W2015\DataAnalysis\Paper3\BottomTrack\';
load([ddir 'BDtrace_floodebb.mat'])
%% Plot Routine
f1 = figure(2);
set(f1,'PaperOrientation','portrait',...
    'position',[400 100   1000  350],...
    'renderer','painters');
names = {'vp1';'vp2';'vp3'};
% symb = {'o';'d';'p'};
lines = {'-';'-';'-'};
tstep = datenum(0,0,0,0,60,0);
hight = datenum(2015,03,09,04,33,00);
sp = zeros(4,1);
cl = flipud([0.1 0.1 0.1;0.5 0.5 0.5;0.7 0.7 0.7]);
pl = zeros(3,1);
space = 300;
% sp(1) = subplot(1,2,1);
for i = 1:3
    time = data.(names{i}).time;
    lvl = data.(names{i}).ble.*1000;
    if i == 1    
        plot(time,zeros(length(time),1),'--k',...
            'linewidth',1.5),hold on
        plot(hight*ones(10,1),linspace(-100,100,10),'-k',...
            'linewidth',1.5)
    end
    plot(time,lvl,lines{i},...
        'color',cl(i,:),'linewidth',1.5)
    lc = 1:space:length(time);
%     pl(i) = plot(time(lc),lvl(lc),...
%         'LineWidth',1.5,...
%         'markeredgecolor','k',...
%         'markerfacecolor',cl(i,:),...
%         'markersize',6);
end
% leg = legend(pl,{'Seaward';'Middle';'Landward'});
dirc = 'd:\Projects\Mekong_W2015\Images\Quadrats\Q2\Reconstructions\';
qdat = 'D:\Projects\Mekong_W2015\Images\Quadrats\Q2\Reconstructions\2A\Analysis\';
fname = 'twobtoc_v2diff.txt';
qname = 'Q2Bstats_fixed.mat';
vpxy = [0.8033 0.245;0.5997 0.245;0.4978 0.245];
rotx = [0,0];
roty = [0,0];
rotz = [0,0];
cc = brewermap(100,'RdYlBu');
colormap(cc);
[z,r] = ImportAsciiRaster([dirc fname]);
Z = (z/10).*1000; %*1000 converts to mm
%     Z = rot90(Z,-1); %rotate Z so it is in the same view as the orthophotos
%create unit vector axis
[m,n] = size(Z);
x = linspace(0,1,n);X = repmat(x,m,1);
y = linspace(1,0,m);Y = repmat(y,n,1)';
%load quadrat data
load([qdat 'Q2Astats_fixed.mat'])

%contour using contourf
sp(2) = subplot(1,2,2);
[c,h] = contourf(X,Y,Z);
set(h,'EdgeColor','none')
hold on
%calculate radii of pneumatophores at level 1 of the reconstruction
nn = length(DATA.layer1.Rcenter);
for j = 1:nn %loop through rows for the Radius and Center points of each layer
    cx = DATA.layer1.Rcenter(j,1);
    cy = DATA.layer1.Rcenter(j,2);
    %transformation matrix, rotate 90deg CCW
    R = [0 -1;1 0];
    rot = R*[cx;cy];
    cx = abs(rot(1));cy = rot(2);
    cr = mean(DATA.layer1.Rradii(j,:));
    ang=0:0.01:2*pi;
    xp=cr.*cos(ang);
    yp=cr.*sin(ang);
    %plot circles with the same radius as the vegetation
    H = patch(cx+xp,cy+yp,1);
    set(H,'FaceColor',[0 0 0],'EdgeColor','none')
end
for j = 1:3
    plot(vpxy(j,1),vpxy(j,2),...
        'Color','k',...
        'MarkerSize',10,...
        'MarkerFaceColor',cl(j,:),...
        'LineWidth',1.5)
end
cb = colorbar;
axis equal
caxis([-30 30])
grid on
set(gca,'gridlinestyle',':')
%Global Plot Adjustments
set(sp(1),...
    'xlim',[time(1) time(end)],...
    'ylim',[-80 40],...
    'ytick',-80:20:40)
datetick(sp(1),'x','HH:MM','keepticks','keeplimits')
set(sp(2),...
    'xlim',[0 1],...
    'ylim',[0 1],...
    'xtick',0:0.2:1,...
    'ytick',0:0.2:1)
set(sp(1),'position',[0.1 0.2 0.45 0.6])
set(sp(2),'position',[0.41 0.18 0.7 0.7])
set(cb,'position',[0.9 0.18 0.025 0.7])
% set(leg,'position',[0.49 0.32 0.02 0.02])
%Labeling
xlabel(sp(1),'Time on 09-03-15 [HH:MM]')
ylabel(sp(1),'Bed Level [mm]')
xlabel(sp(2),'Across-Shore Distance [m]')
ylabel(sp(2),'Along-Shore Distance [m]')
ylabel(cb,'Bed Level Change [mm]')
prettyfigures('text',13,'labels',14,'box',1)
sfdir = 'g:\GradSchool\DataAnalysis\Paper3\Figures\';
% export_fig([sfdir 'BDtrace_bdmap_v2'],'-pdf','-nocrop')