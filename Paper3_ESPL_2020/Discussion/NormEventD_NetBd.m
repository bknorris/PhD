%Plot the combined event data
clear,close all
load('f:\Mekong_W2015\DataAnalysis\Paper3\Wavelet\CombEventdata_flood.mat');
data.flood = dat;
load('f:\Mekong_W2015\DataAnalysis\Paper3\Wavelet\CombEventdata_ebb.mat');
data.ebb = dat;clear dat
fn = fieldnames(data);
sfdir = 'g:\GradSchool\DataAnalysis\Paper3\Figures\';
flds = {'mud';'fringe';'forest'};
band = {'wave';'ig'};
cl = [178,24,43;5,55,122]./255;
symb = {'o';'d';};
%% Net Bed Level - Event Length (Acc & Ero)
sp = zeros(3,2);
eb = zeros(2,1);
f1 = figure(1);
set(f1,'PaperOrientation','portrait',...
    'position',[400 100   800   500],...
    'renderer','painters');
w = [1 4; 2 5; 3 6];
tsep = [-0.001;0.001];
for i = 1:2
    disp(fn{i})
    for j = 1:length(flds)
        disp(flds{j})
        %Accretion
        disp('Accretion')
        sp(w(j,1)) = subplot(2,3,w(j,1));
        deltbd = [data.(fn{i}).(flds{j}).wave.deltbd; data.(fn{i}).(flds{j}).ig.deltbd];
        id = find(deltbd>0);
        eventl = [data.(fn{i}).(flds{j}).wave.eventl; data.(fn{i}).(flds{j}).ig.eventl];
        ntime = eventl./max(eventl);
        xs = ntime(id);ys = deltbd(id);
        xs(ys>0.1) = [];ys(ys>0.1) = [];
        pf = polyfit(xs,ys,1);
        xx = 0:0.1:1;
        pv = polyval(pf,xx);
        %rsquared
        rsq = corrcoef(xs,ys);rsq = rsq(1,2);
        plot(xs,ys,symb{i},'color',cl(i,:),'markerfacecolor',cl(i,:),'markersize',4),hold on
%         plot(xx,pv,'-',...
%             'linewidth',2.5,...
%             'color','w')
%         plot(xx,pv,'-',...
%             'linewidth',1.5,...
%             'color',cl(i,:))
%         rtext = sprintf('R^2 = %0.2f',rsq);
%         text(xx(5),pv(5)+tsep(i),rtext,'color',cl(i,:))
        fprintf('Median: %0.2f\n',median(deltbd(deltbd>0.01|deltbd<-0.01))*1000)
        fprintf('Max: %0.2f\n',max(ys)*1000)
        fprintf('25th and 75th (eventl) quantiles: %0.2f and %0.2f seconds\n',quantile(xs.*max(eventl)*60,[0.25 0.75]))
        fprintf('25th and 75th (eventl - IG only) quantiles: %0.2f and %0.2f min\n',quantile(data.(fn{i}).(flds{j}).ig.eventl,[0.25 0.75]))
        fprintf('25th and 75th (deltbd) quantiles: %0.2d and %0.2d mm\n',quantile(ys.*1000,[0.25 0.75]))
        %Erosion
        disp('Erosion')
        sp(w(j,2)) = subplot(2,3,w(j,2));
        id = find(deltbd<0);
        xs = ntime(id);ys = deltbd(id);
        xs(ys<-0.1) = [];ys(ys<-0.1) = [];
        pf = polyfit(xs,ys,1);
        xx = 0:0.1:1;
        pv = polyval(pf,xx);
        rsq = corrcoef(xs,ys);rsq = rsq(1,2);
        eb(i) = plot(xs,ys,symb{i},'color',cl(i,:),'markerfacecolor',cl(i,:),'markersize',4);hold on
%         plot(xx,pv,'-',...
%             'linewidth',2.5,...
%             'color','w')
%         plot(xx,pv,'-',...
%             'linewidth',1.5,...
%             'color',cl(i,:))
%         rtext = sprintf('R^2 = %0.2f',rsq);
%         text(xx(5),pv(5)+tsep(i),rtext,'color',cl(i,:))
        fprintf('Min: %0.2f\n',min(ys)*1000)
        fprintf('25th and 75th (eventl) quantiles: %0.2f and %0.2f seconds\n',quantile(xs.*max(eventl)*60,[0.25 0.75]))
        fprintf('25th and 75th (eventl - IG only) quantiles: %0.2f and %0.2f min\n',quantile(data.(fn{i}).(flds{j}).ig.eventl,[0.25 0.75]))
        fprintf('25th and 75th (deltbd) quantiles: %0.2d and %0.2d mm\n',quantile(ys*1000,[0.25 0.75]))    
    end
end
%Plot adjustments
set(sp,'xlim',[-0.02 1.05],'xtick',0:0.2:1)
%Top row
set(sp(1),'ylim',[-0.005 0.1],'position',[0.12 0.55 0.25 0.37])
set(sp(2),'ylim',[-0.005 0.1],'position',[0.405 0.55 0.25 0.37])
set(sp(3),'ylim',[-0.005 0.1],'position',[0.69 0.55 0.25 0.37])
title(sp(1),'Mudflat')
title(sp(2),'Fringe')
title(sp(3),'Forest')
suplabel('Net elevation change over event [m]','y');
%Bottom row
set(sp(4),'ylim',[-0.1 0.005],'position',[0.12 0.14 0.25 0.37])
set(sp(5),'ylim',[-0.1 0.005],'position',[0.405 0.14 0.25 0.37])
set(sp(6),'ylim',[-0.1 0.005],'position',[0.69 0.14 0.25 0.37])
xlabel(sp(4),'Normalized Time')
xlabel(sp(5),'Normalized Time')
xlabel(sp(6),'Normalized Time')
set([sp(2) sp(3) sp(5) sp(6)],'yticklabel',[])
set([sp(1) sp(2) sp(3)],'xticklabel',[])
leg = legend([eb(1),eb(2)],{'Flood Tide';'Ebb Tide'});
set(leg,'position',[0.77 0.19 0.005 0.005])
prettyfigures('text',12,'labels',13,'box',1)
% export_fig([sfdir 'EventlNetChange_v2'],'-pdf','-nocrop')
 
