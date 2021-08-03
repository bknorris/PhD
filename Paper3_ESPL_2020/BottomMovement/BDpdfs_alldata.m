%plot the pdf of the bottom trace relative to the initial bed elevation
clear,close all
load('f:\Mekong_W2015\DataAnalysis\Paper3\Wavelet\CombEventData.mat');
fn = fieldnames(dat);
f1 = figure(1);
set(f1,'PaperOrientation','portrait',...
    'position',[400 100   800   500],...
    'renderer','painters');
cl = [0.7 0.7 0.7;0.4 0.4 0.4;0 0 0];
symb = {'o','.','^'};
p = zeros(1,3);
sp(1) = subplot(211);
for i = 1:3
    bdn = dat.(fn{i}).bdmed;
    bins = linspace(-0.1,0.1,20);
    N = hist(bdn,bins);
    p(i) = plot(bins,(N./numel(bdn))*100,...
        '-',...
        'marker',symb{i},...
        'color',cl(i,:),...
        'linewidth',1.5);hold on
end
sp(2) = subplot(212);
for i = 1:3
    deltabd = dat.(fn{i}).deltbd;
    bins = linspace(-0.1,0.1,20);
    N = hist(deltabd,bins);
    p(i) = plot(bins,(N./numel(deltabd))*100,...
        '-',...
        'marker',symb{i},...
        'color',cl(i,:),...
        'linewidth',1.5);hold on
end
set(sp(1),'position',[0.13 0.58 0.8 0.35],...
    'ylim',[-2 50],'xlim',[-0.06 0.06])
set(sp(2),'position',[0.13 0.11 0.8 0.35],...
    'ylim',[-2 50],'xlim',[-0.04 0.04])
xlabel(sp(1),'Elevation relative to initial (m)')
xlabel(sp(2),'Net elevation change over event (m)')
ylabel(sp(1),'Percent occurence')
ylabel(sp(2),'Percent occurence')
leg = legend(p,{'Mudflat','Fringe','Forest'});
prettyfigures('text',12,'labels',13,'box',1)
set(leg,'position',[0.83 0.35 0.05 0.05],'box','off')
sfdir = 'd:\Projects\Mekong_W2015\Figures\Paper3\Discussion\';
% export_fig([sfdir 'AllData_bdevel_pdf'],'-png')