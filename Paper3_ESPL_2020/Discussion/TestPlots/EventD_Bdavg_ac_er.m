clear, close all
load('g:\Mekong_W2015\DataAnalysis\Paper3\Wavelet\12-03-15\waveIGevents.mat')
sfdir = 'f:\GradSchool\DataAnalysis\Paper3\WorkingFigures\Wavelet\';
expname = 'F2F3_2';
exploc = {'Mudflat';'Fringe';'Forest'};
% exploc = {'x = -10cm';'x = 10cm';'x = 20cm'};

ltxt = cell(3,1);
p1 = zeros(3,1);
fn = fieldnames(eventd);
cl = [207 176 126;
    60 166 74;
    4 76 41]./255;
symb = {'o';'^';'s'};
f1 = figure(1);
set(f1,'PaperOrientation','portrait',...
    'position',[400 100   800   400],...
    'renderer','painters')
for i = 1:length(fn);
    %normalize event length by longest event, bin median bed level by event
    %length
    eventl = [eventd.(fn{i}).wave.eventl; eventd.(fn{i}).ig.eventl];
    eventl = eventl(~isnan(eventl));
    nevtl = linspace(min(eventl./max(eventl)),max(eventl./max(eventl)),30);
    bdmed = [eventd.(fn{i}).wave.deltbd; eventd.(fn{i}).ig.deltbd];
    bdmed = bdmed(~isnan(bdmed));
    %%%
    for ii = 1:2
        sp(ii) = subplot(1,2,ii); %positive bd
        if ii == 1
            id = bdmed>0;
        else
            id = bdmed<0;
        end
        [b,n,s] = bindata(eventl(id)./max(eventl(id)),bdmed(id),nevtl);
        p1(i) = errorbar(nevtl,b,s,symb{i},...
            'color',cl(i,:),...
            'linewidth',1.5,...
            'markersize',8,...
            'markeredgecolor','k',...
            'markerfacecolor',cl(i,:));hold on
    end
    ltxt{i} = exploc{i};
end
leg = legend(p1(p1~=0),ltxt(~cellfun('isempty',ltxt)));
%Global plot adjustments
set(sp(1),'ylim',[-1E-4 60E-3],'xlim',[-0.05 1.05],...
    'position',[0.12 0.15 0.35 0.75])
set(sp(2),'ylim',[-60E-3 1E-4],'xlim',[-0.05 1.05],...
    'position',[0.57 0.15 0.35 0.75])
ylabel(sp(1),'Net Bed Level Change (m)')
xlabel(sp(1),'Normalized Event Length')
xlabel(sp(2),'Normalized Event Length')
title(sp(1),'Accretion')
title(sp(2),'Erosion')
prettyfigures('text',12,'labels',13,'box',1)
set(leg,'position',[0.61,0.17,0.1 0.2])
export_fig([sfdir expname '_evlBDavg_acc_ero'],'-png')

