%Plot revised versions of the bed level-forcing variable plots
%15/04/20
%
clear,close all
%define working paths
dpath = '\DataAnalysis\Paper3\';
ypath = 'Mekong_W2015';
%load run file to tell program which files & VPs to load
rdir = 'e:\MemoryStick\GradSchool\DataAnalysis\Paper3\ExperimentalDesign\';
fid = fopen([rdir 'AutoForcingPlotrunfile.csv']);
rfile = textscan(fid,'%s%s%s%s%n%n%n','delimiter',',');
rdate = rfile{1};
rvp = rfile{2};
rtide = rfile{3};
rexp = rfile{4};
rvph = rfile{5};
ra = rfile{6};
rd = rfile{7};

%Set up plot parameters
symb = {'^';'o';'s';'d'};
color = brewermap(4,'greys');
tauc = [0.26 0.32]; %critical shear stress

%Load the data
npath = ['e:\' ypath dpath];
folders = dir([npath '\VPs\']);
folders = {folders(3:end).name};
c = 1;
for ii = 1:length(folders)
    if any(strcmp(folders{ii},rdate))
        vpid = strcmp(folders{ii},rdate);vpid = find(vpid);
        disp(['Loading files from ' npath 'Wavelet\' folders{ii} '\'])
        load([npath 'Wavelet\' folders{ii} '\','waveIGevents_v4.mat']);
    else
        continue
    end
    f1 = figure(1);
    set(f1,'PaperOrientation','portrait',...
        'position',[400 100   850   450]);
    sp(1) = subplot(121);
    taub = eventd.(rvp{vpid(1)}).ig.taub;
    eps = eventd.(rvp{vpid(1)}).ig.eps;
    bdmean = (rvph(vpid(1))-eventd.(rvp{vpid(1)}).ig.bdmean);
    idx = find(taub>tauc(1));
    p(c) = plot(decimate(eps(idx),10),decimate(bdmean(idx),10),...
        symb{c},'color','k',...
        'markerfacecolor',color(c,:),...
        'linewidth',1.5,...
        'markersize',6); hold on
    
    sp(2) = subplot(122);
    taub = eventd.(rvp{vpid(2)}).ig.taub;
    eps = abs(eventd.(rvp{vpid(2)}).ig.eps);
    bdmean = (rvph(vpid(2))-eventd.(rvp{vpid(2)}).ig.bdmean);
    idx = find(taub>tauc(2));
    p(c) = plot(decimate(eps(idx),8),decimate(bdmean(idx),8),...
        symb{c},...
        'color','k',...
        'markerfacecolor',color(c,:),...
        'linewidth',1.5,...
        'markersize',6); hold on
    c = c+1;
end
legtext = {'Ebb tide, n = 102';'Flood tide, n = 102';'Flood tide, n = 45';'Ebb tide, n = 45'};
leg1 = legend(sp(1),legtext);
legtext = {'Ebb tide, n = 84';'Flood tide, n = 84';'Flood tide, n = 37';'Ebb tide, n = 37'};
leg2 = legend(sp(2),legtext);
set(sp(1),'xlim',[-3E-6 2.5E-4],...
    'yticklabel',[-20 -15 -10 -5 0 5 10 15 20 25],...
    'position',[0.12 0.12 0.38 0.8])
ytick = get(sp(1),'ytick');
set(sp(2),'xlim',[-1E-6 8E-5],...
    'ytick',ytick,...
    'yticklabel',[],...
    'position',[0.56 0.12 0.38 0.8])
xlabel(sp(1),'\epsilon [m^2/s^3]')
xlabel(sp(2),'\epsilon [m^2/s^3]')
ylabel(sp(1),'Mean BLE During Event [mm]')
title(sp(1),'Fringe Experiments')
title(sp(2),'Forest Experiments')
prettyfigures('text',12,'labels',13,'box',1,'gcolor','k')
set(f1,'units','inches');
pos = get(f1,'Position');
set(f1,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
sfdir = 'e:\MemoryStick\GradSchool\DataAnalysis\Paper3\Figures\';
export_fig([sfdir 'TurbDissip_BLE_Fringe_forest'],'-pdf','-nocrop')

