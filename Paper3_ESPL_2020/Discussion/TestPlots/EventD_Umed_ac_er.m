clear, close all
load('D:\Mekong_F2014\DataAnalysis\Paper3\Wavelet\30-09-14\waveIGevents.mat')
sfdir = 'e:\GradSchool\DataAnalysis\Paper3\WorkingFigures\Wavelet\';
expname = 'F2F_2';
exploc = {'Mudflat';'Fringe';'Forest'};
ltxt = cell(3,1);
p1 = zeros(3,1);
p2 = zeros(3,1);
fn = fieldnames(eventd);
cl = [207 176 126;
    60 166 74;
    4 76 41]./255;
cc = [180 180 180;124 124 124;51 51 51]./255;
symb = {'o';'^';'s'};
f1 = figure(1);
    set(f1,'PaperOrientation','portrait',...
        'position',[400 100   800   400],...
        'renderer','painters')
for i = 1:length(fn)
    %normalize event length by longest event, bin median velocity by event
    %length
    deltbd = [eventd.(fn{i}).wave.deltbd; eventd.(fn{i}).ig.deltbd];
    deltbd = deltbd(~isnan(deltbd));
    eventl = [eventd.(fn{i}).wave.eventl; eventd.(fn{i}).ig.eventl];
    eventl = eventl(~isnan(eventl));
    nevtl = linspace(min(eventl./max(eventl)),max(eventl./max(eventl)),30);
    umed = [eventd.(fn{i}).wave.umed; eventd.(fn{i}).ig.umed]; 
    umed = umed(~isnan(umed));
    %%%
    sp(1) = subplot(121); %positive bd
    id = deltbd>0;
    [b,~,s] = bindata(eventl(id)./max(eventl(id)),umed(id),nevtl);
    errorbar(nevtl,b,s,symb{i},...
        'color',cl(i,:),...
        'linewidth',1.5,...
        'markersize',8,...
        'markeredgecolor','k',...
        'markerfacecolor',cl(i,:)),hold on
    %%%
    sp(2) = subplot(122); %positive bd
    id = deltbd<0;
    [b,~,s] = bindata(eventl(id)./max(eventl(id)),umed(id),nevtl);
    p1(i) = errorbar(nevtl,b,s,symb{i},...
        'color',cl(i,:),...
        'linewidth',1.5,...
        'markersize',8,...
        'markeredgecolor','k',...
        'markerfacecolor',cl(i,:));hold on
    ltxt{i} = exploc{i};
end
leg = legend(p1(p1~=0),ltxt(~cellfun('isempty',ltxt)));
%Global plot adjustments
set(sp(1),'ylim',[-6E-2 6E-2],'xlim',[-0.05 1.05],...
    'position',[0.15 0.15 0.35 0.75])
set(sp(2),'ylim',[-6E-2 6E-2],'xlim',[-0.05 1.05],...
    'yticklabel',[],...
    'position',[0.55 0.15 0.35 0.75])
ylabel(sp(1),'Median velocity (m/s)')
xlabel(sp(1),'Normalized Event Length')
xlabel(sp(2),'Normalized Event Length')
title(sp(1),'Accretion')
title(sp(2),'Erosion')
prettyfigures('text',12,'labels',13,'box',1)
export_fig([sfdir expname '_evlUmed_acc_ero'],'-png','-nocrop')

