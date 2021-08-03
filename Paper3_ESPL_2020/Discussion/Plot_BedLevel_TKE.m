%Plot revised versions of the bed level-forcing variable plots
%15/04/20
%
clear,close all
%define working paths
dpath = '\DataAnalysis\Paper3\';
ypath = 'Mekong_W2015';
%load run file to tell program which files & VPs to load
rdir = 'f:\MemoryStick\GradSchool\DataAnalysis\Paper3\ExperimentalDesign\';
fid = fopen([rdir 'AutoForcingPlotrunfile_v2.csv']);
rfile = textscan(fid,'%s%s%s%s%n%n%n','delimiter',',');
rdate = rfile{1};
rvp = rfile{2};
rtide = rfile{3};
rexp = rfile{4};
rvph = rfile{5};
ra = rfile{6};
rd = rfile{7};

%Set up plot parameters
symb = {'^';'o';'d'};
c = [207 176 126;
    60 166 74;
    4 76 41]./255;
tauc = [0.18 0.26 0.32]; %critical shear stress

%Load the data
npath = ['f:\' ypath dpath];
folders = dir([npath '\VPs\']);
folders = {folders(3:end).name};
for ii = 1:length(folders)
    if any(strcmp(folders{ii},rdate))
        vpid = strcmp(folders{ii},rdate);vpid = find(vpid);
        disp(['Loading files from ' npath 'Wavelet\' folders{ii} '\'])
        load([npath 'Wavelet\' folders{ii} '\','waveIGevents_v5.mat']);
    else
        continue
    end
    f1 = figure(1);
    set(f1,'PaperOrientation','portrait',...
        'position',[400 100   750   350]);
    sp(1) = subplot(121);
    plot(linspace(0.001,0.1,10),zeros(1,10),':k','linewidth',1.5), hold on
    for i = 1:3
    u = abs(eventd.(rvp{vpid(i)}).ig.usqd);
    bdmean = eventd.(rvp{vpid(i)}).ig.deltbd.*1000;
    p(i) = plot(u,bdmean,...
        symb{i},'color','k',...
        'markerfacecolor',c(i,:),...
        'linewidth',1.5,...
        'markersize',7);
    end 
    sp(2) = subplot(122);
    plot(linspace(0.5E-5,1E-3,10),zeros(1,10),':k','linewidth',1.5), hold on
    for i = 1:3
    eps = eventd.(rvp{vpid(i)}).ig.eps20;
    bdmean = eventd.(rvp{vpid(i)}).ig.deltbd.*1000;
    p(i) = plot(eps,bdmean,...
        symb{i},'color','k',...
        'markerfacecolor',c(i,:),...
        'linewidth',1.5,...
        'markersize',7);
    end
end
legtext = {'Mudflat';'Fringe';'Forest'};
leg1 = legend(p,legtext,'box','off');
% legtext = {'Ebb tide, n = 84';'Flood tide, n = 84';'Flood tide, n = 37';'Ebb tide, n = 37'};
% leg2 = legend(sp(2),legtext);
set(sp(1),'xlim',[0.001 0.1],...
    'ylim',[-3 4],...
    'position',[0.15 0.18 0.35 0.76])
set(sp(2),'xlim',[5E-5 1E-3],...
    'ylim',[-3 4],...
    'ytick',[],...
    'position',[0.55 0.18 0.35 0.76])
xlabel(sp(1),'$\bar{u} {[m/s]}$','interpreter','latex')
xlabel(sp(2),'k [m^2/s^2]')
% xlabel(sp(2),'\epsilon [m^2/s^3]')
ylabel(sp(1),'Net BLE During Event [mm]')
% title(sp(2),'Forest Experiments')
prettyfigures('text',12,'labels',13,'box',1,'gcolor','k')
set(f1,'units','inches');
% pos = get(f1,'Position');
% set(f1,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
% sfdir = 'd:\MemoryStick\GradSchool\DataAnalysis\Paper3\Figures\';
% export_fig([sfdir 'U_TKE_BLE_v2'],'-pdf','-nocrop')
% 
