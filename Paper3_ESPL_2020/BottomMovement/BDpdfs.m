%plot the pdf of the bottom trace relative to the initial bed elevation
clear,close all
bd = load('d:\Projects\Mekong_W2015\DataAnalysis\Paper3\BottomTrack\F2F3_1_bdtrace.mat');
load('G:\Mekong_W2015\DataAnalysis\Paper3\Wavelet\11-03-15\waveIGevents.mat')
fn = fieldnames(bd);
f1 = figure(1);
set(f1,'PaperOrientation','portrait',...
    'position',[400 100   800   500],...
    'renderer','painters');
cl = [0.7 0.7 0.7;0.4 0.4 0.4;0 0 0];
line = {'-';'--';':'};
p = zeros(1,3);
sp(1) = subplot(211);
for i = 1:3
    bdist = bd.(fn{i}).bdist;
    bdist = bdist(~isnan(bdist));
    bint = bdist(1);
    bdn = bdist-bint;
    bins = linspace(min(bdn),max(bdn),20);
    N = hist(bdn,bins);
    p(i) = plot(bins,(N./numel(bdn))*100,...
        line{i},...
        'color',cl(i,:),...
        'linewidth',1.5);hold on
end
sp(2) = subplot(212);
for i = 1:3
    deltabd = [eventd.(fn{i}).wave.deltbd; eventd.(fn{i}).ig.deltbd];
    deltabd = deltabd(~isnan(deltabd));
    %deltabd has some huge outliers, remove
    id = find(deltabd > nanmean(deltabd)-2*nanstd(deltabd)...
        & deltabd < nanmean(deltabd)+2*nanstd(deltabd));
    bins = linspace(min(deltabd(id)),max(deltabd(id)),20);
    N = hist(deltabd(id),bins);
    p(i) = plot(bins,(N./numel(deltabd(id)))*100,...
        line{i},...
        'color',cl(i,:),...
        'linewidth',1.5);hold on
end
set(sp(1),'position',[0.13 0.58 0.8 0.35],...
    'ylim',[0 40])
set(sp(2),'position',[0.13 0.11 0.8 0.35],...
    'ylim',[0 100],'xlim',[-2E-3 2E-3])
xlabel(sp(1),'Elevation relative to initial (m)')
xlabel(sp(2),'Net elevation change over event (m)')
ylabel(sp(1),'Percent occurence')
ylabel(sp(2),'Percent occurence')
leg = legend(p,{'Mudflat','Fringe','Forest'});
prettyfigures('text',12,'labels',13,'box',1)
set(leg,'position',[0.83 0.35 0.05 0.05],'box','off')
sfdir = 'd:\Projects\Mekong_W2015\Figures\Paper3\Discussion\';
% export_fig([sfdir 'F2F3_bdevel_pdf'],'-png')