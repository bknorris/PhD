%Plot the combined event data
clear,close all
load('e:\Mekong_W2015\DataAnalysis\Paper3\Wavelet\CombEventdata_flood_v4.mat');
data.flood = dat;
load('e:\Mekong_W2015\DataAnalysis\Paper3\Wavelet\CombEventdata_ebb_v4.mat');
data.ebb = dat;clear dat
fn = fieldnames(data);
sfdir = 'd:\GradSchool\DataAnalysis\Paper3\Figures\';
flds = {'mud';'fringe';'forest'};
band = {'wave';'ig'};
cl = [0.8 0.8 0.8; 0.2 0.2 0.2];
symb = {'o';'^';'s'};
%% Net Bed Level - Event Length (Acc & Ero)
sp = zeros(3,2);
eb = zeros(3,2);
f1 = figure(1);
set(f1,'PaperOrientation','portrait',...
    'position',[400 100   800   500],...
    'renderer','painters');
w = [1 4; 2 5; 3 6];
for i = 1:2
    for j = 1:length(flds)
        %Accretion
        sp(w(j,1)) = subplot(2,3,w(j,1));
        deltbd = [data.(fn{i}).(flds{j}).wave.deltbd; data.(fn{i}).(flds{j}).ig.deltbd];
        id = find(deltbd>0);
        eventl = [data.(fn{i}).(flds{j}).wave.eventl; data.(fn{i}).(flds{j}).ig.eventl];
        ntime = eventl./max(eventl);
        bins = linspace(min(ntime(id)),max(ntime(id)),13);
        [b,~,q1,q3] = binmedian(ntime(id),deltbd(id),bins);
        eb(w(j,1)) = errorbar(bins,b,q1,q3,symb{j},...
            'color','k',...
            'linewidth',1.5,...
            'capsize',3,...
            'markersize',5,...
            'markeredgecolor','k',...
            'markerfacecolor',cl(i,:));hold on
        %Erosion
        sp(w(j,2)) = subplot(2,3,w(j,2));
        id = find(deltbd<0);
        bins = linspace(min(ntime(id)),max(ntime(id)),13);
        [b,~,q1,q3] = binmedian(ntime(id),deltbd(id),bins);
        eb(w(j,2)) = errorbar(bins,b,q1,q3,symb{j},...
            'color','k',...
            'linewidth',1.5,...
            'capsize',3,...
            'markersize',5,...
            'markeredgecolor','k',...
            'markerfacecolor',cl(i,:));hold on
    end
end
leg = legend([eb(3,1),eb(3,2)],{'Flood Tide';'Ebb Tide'});
set(leg,'position',[0.92 0.89 0.005 0.005])
%Plot adjustments
set(sp,'xlim',[-0.02 1.05],'xtick',0:0.2:1)
%Top row
set(sp(1),'ylim',[-0.005 0.1],'position',[0.12 0.59 0.23 0.35])
set(sp(2),'ylim',[-0.004 0.08],'position',[0.42 0.59 0.23 0.35])
set(sp(3),'ylim',[-0.002 0.04],'position',[0.72 0.59 0.23 0.35])
title(sp(1),'Mudflat')
title(sp(2),'Fringe')
title(sp(3),'Forest')
suplabel('Net elevation change over event [m]','y');
%Bottom row
set(sp(4),'ylim',[-0.1 0.005],'position',[0.12 0.14 0.23 0.35])
set(sp(5),'ylim',[-0.08 0.004],'position',[0.42 0.14 0.23 0.35])
set(sp(6),'ylim',[-0.04 0.002],'position',[0.72 0.14 0.23 0.35])
xlabel(sp(4),'Normalized Time')
xlabel(sp(5),'Normalized Time')
xlabel(sp(6),'Normalized Time')
prettyfigures('text',12,'labels',13,'box',1)
% export_fig([sfdir 'NetChangeEventl'],'-pdf','-nocrop')
%% Net Bed Level - Epsilon (Acc & Ero)
sp = zeros(5,2);
eb = zeros(3,1);
f2 = figure(2);
set(f2,'PaperOrientation','portrait',...
    'position',[400 100   1000   500],...
    'renderer','painters');
w = [1 6; 2 7; 3 8;4 9;5 10];
cl = [207 176 126;
    60 166 74;
    4 76 41]./255;
tcr = [0.32,0.4,1.13]; %from Critical_sed_params.m
lspec = {'-';'--';':'};
fields = {'depth';'sigh';'usqd';'taub';'eps'};
for i = 1:length(fields)
    for j = 1:length(flds)
        %Flood Tide
        sp(w(i,1)) = subplot(2,5,w(i,1));
        bdmed = [data.flood.(flds{j}).wave.bdmed; data.flood.(flds{j}).ig.bdmed];
        xx = [data.flood.(flds{j}).wave.(fields{i}); data.flood.(flds{j}).ig.(fields{i})];
        bins = linspace(min(xx),max(xx),15);
        [b,~,q1,q3] = binmedian(xx,bdmed,bins);
        %interp nans
        q1(q1==0)=NaN;q3(q3==0)=NaN;
        b = interpnan(b);q1 = interpnan(q1);q3 = interpnan(q3);
        if j == 1
            plot(linspace(0,100,10),zeros(1,10),'-k','linewidth',1.5);hold on
        end
        bl = boundedline(bins,b,[abs(q1) abs(q3)],'alpha','transparency',0.3,...
            'cmap',cl(j,:));hold on
        set(bl,'linewidth',1.5)
        if i == 1,eb(j) = bl;end
        if strcmp(fields{i},'taub')
            for k = 1:3
            plot(ones(10,1)*tcr(k),linspace(-0.1,0.1,10),lspec{k},...
                'linewidth',1,'color','k')
            end
        end
        %Ebb Tide
        sp(w(i,2)) = subplot(2,5,w(i,2));
        bdmed = [data.ebb.(flds{j}).wave.bdmed; data.ebb.(flds{j}).ig.bdmed];
        xx = [data.ebb.(flds{j}).wave.(fields{i}); data.ebb.(flds{j}).ig.(fields{i})];
        bins = linspace(min(xx),max(xx),15);
        [b,~,q1,q3] = binmedian(xx,bdmed,bins);
        %interp nans
        q1(q1==0)=NaN;q3(q3==0)=NaN;
        b = interpnan(b);q1 = interpnan(q1);q3 = interpnan(q3);
        if j == 1
            plot(linspace(0,100,10),zeros(1,10),'-k','linewidth',1.5);hold on
        end
        bl = boundedline(bins,b,[abs(q1) abs(q3)],'alpha','transparency',0.3,...
            'cmap',cl(j,:));hold on
        set(bl,'linewidth',1.5)
        if strcmp(fields{i},'taub')
            for k = 1:3
                plot(ones(10,1)*tcr(k),linspace(-0.1,0.1,10),lspec{k},...
                    'linewidth',1,'color','k')
            end
        end
    end
end
leg = legend(eb,{'Mudflat';'Fringe';'Forest'});
set(leg,'position',[0.94 0.87 0.005 0.005])
% Plot adjustments
% Top row
set([sp(2) sp(3) sp(4) sp(5) sp(7) sp(8) sp(9) sp(10)],'yticklabel',[])
set([sp(1) sp(2) sp(3) sp(4) sp(5)],'xticklabel',[])
set(sp(1),'xlim',[0 1.5],'ylim',[-0.1 0.1],'position',[0.11 0.59 0.15 0.35])
set(sp(2),'xlim',[0 0.6],'xtick',0:0.2:0.6,'ylim',[-0.1 0.1],'position',[0.29 0.59 0.15 0.35])
set(sp(3),'xlim',[0 0.3],'xtick',0:0.1:0.3,'ylim',[-0.1 0.1],'position',[0.47 0.59 0.15 0.35])
set(sp(4),'xlim',[0 1.5],'ylim',[-0.1 0.1],'position',[0.65 0.59 0.15 0.35])
set(sp(5),'xlim',[0 1E-3],'ylim',[-0.1 0.1],'position',[0.83 0.59 0.15 0.35])
suplabel('Median elev. relative to event beginning [m]','y');
%Bottom row
set(sp(6),'xlim',[0 1.5],'ylim',[-0.1 0.1],'position',[0.11 0.15 0.15 0.35])
set(sp(7),'xlim',[0 0.6],'xtick',0:0.2:0.6,'ylim',[-0.1 0.1],'position',[0.29 0.15 0.15 0.35])
set(sp(8),'xlim',[0 0.3],'xtick',0:0.1:0.3,'ylim',[-0.1 0.1],'position',[0.47 0.15 0.15 0.35])
set(sp(9),'xlim',[0 1.5],'ylim',[-0.1 0.1],'position',[0.65 0.15 0.15 0.35])
set(sp(10),'xlim',[0 1E-3],'ylim',[-0.1 0.1],'position',[0.83 0.15 0.15 0.35])
xlabel(sp(6),'Water Depth [m]')
xlabel(sp(7),'H_s [m]')
xlabel(sp(8),'u^2 [m^2s^{-2}]')
xlabel(sp(9),'\tau_b [Pa]')
xlabel(sp(10),'\epsilon [Wkg^{-1}]')
prettyfigures('text',12,'labels',13,'box',1)
set(f2,'units','inches');
pos = get(f2,'Position');
set(f2,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
% print(f2,[sfdir 'ForcingsBdmed'],'-dpdf','-r0')
