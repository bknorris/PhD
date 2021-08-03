%plot a figure showing a rose of u,v to get the 'sector' of current
%directions to plot on the wind rose
clear, close all
aqfil = {'F2F3_11_AD5117.mat';'FSS_08_AD5116.mat'};
aqdir = {'e:\Mekong_W2015\DataAnalysis\Paper3\WaveStats\11-03-15\';...
    'e:\Mekong_W2015\DataAnalysis\Paper3\WaveStats\08-03-15\'};
sdir = 'd:\GradSchool\DataAnalysis\Paper3\Figures\';
%Load the VP positions data
f1 = figure(1);
set(f1,'PaperOrientation','portrait',...
    'position',[400 200   800   550]);
set(gcf,'color','w','paperpositionmode','auto',...
    'renderer','painters')
for k = 1:2
    %% Load the Data
    disp(['Loading ' aqdir{k} aqfil{k}])
    load([aqdir{k} aqfil{k}])
    %average u,v vels
    u = nanmean(wave.u,2);
    v = nanmean(wave.v,2);
    um = runningmean(u,480);
    vm = runningmean(v,480);

    
    %Plot AQDP wave rose
    sp(k) = subplot(2,1,k);
    cmap = brewermap(100,'PuBu');
    [spd,dir] = cmguv2spd(um,vm);
    h = WindRose(dir,spd,...
        'anglenorth',0,...
        'angleeast',90,...
        'titlestring',[],...
        'legendtype',1,...
        'axes',sp(k));
end
