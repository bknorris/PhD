%Plot bed level depth dependency - GM 2018 Norris et al.
clear,close all
load('d:\Projects\Mekong_W2015\DataAnalysis\Paper3\Wavelet\CombEventdata_flood.mat');
data.flood = dat;
load('d:\Projects\Mekong_W2015\DataAnalysis\Paper3\Wavelet\CombEventdata_ebb.mat');
data.ebb = dat;clear dat
fn = fieldnames(data);
sfdir = 'g:\GradSchool\DataAnalysis\Paper3\Figures\';
flds = {'mud';'fringe';'forest'};
f1 = figure(1);
set(f1,'PaperOrientation','portrait',...
    'position',[400 100   700   400],...
    'renderer','painters');
symb = {'o';'^';'s'};
cl = [207 176 126;
    60 166 74;
    4 76 41]./255;
sp = zeros(2,1);
eb = zeros(3,1);
for i = 1:2
    for j = 1:length(flds)
        sp(i) = subplot(1,2,i);
        bdmed = [data.(fn{i}).(flds{j}).wave.bdmed; data.(fn{i}).(flds{j}).ig.bdmed];
        H = [data.(fn{i}).(flds{j}).wave.depth; data.(fn{i}).(flds{j}).ig.depth];
        if j == 3
            H(bdmed<-0.06) = [];bdmed(bdmed<-0.06) = [];
        else
            H(bdmed<-0.1) = [];bdmed(bdmed<-0.1) = [];
        end
        id = find(bdmed<-0.003 | bdmed > 0.003);
        bins = linspace(min(H(id)),max(H(id)),12);
        [b,~,q1,q3] = binmedian(H(id),bdmed(id),bins);
        %interp nans
        q1(q1==0)=NaN;q3(q3==0)=NaN;
        b = interpnan(b);q1 = interpnan(q1);q3 = interpnan(q3);
        if j == 1
            plot(linspace(0,2,10),zeros(1,10),'-k','linewidth',1.5);hold on
        end
        bl = boundedline(bins,b,[abs(q1) abs(q3)],'alpha','transparency',0.4,...
            'cmap',cl(j,:));hold on
        set(bl,'linewidth',1.5)
        if i == 2
            eb(j) = bl;
        end
% plot(H,bdmed,'.','color',cl(j,:)),hold on
    end
end
%Calculate r-squared values
rsq = zeros(3,2);
line = {'-';'--';'-.'};
for i = 1:2
    disp([fn{i} ' tidal stage'])
    for j = 1:length(flds)
        disp(flds{j})
        bdmed = [data.(fn{i}).(flds{j}).wave.bdmed; data.(fn{i}).(flds{j}).ig.bdmed];
        H = [data.(fn{i}).(flds{j}).wave.depth; data.(fn{i}).(flds{j}).ig.depth];
        if j == 3
            H(bdmed<-0.06) = [];bdmed(bdmed<-0.06) = [];
        else
            H(bdmed<-0.1) = [];bdmed(bdmed<-0.1) = [];
        end       
        id = find(bdmed<-0.003 | bdmed > 0.003);
        pf = polyfit(H(id),bdmed(id),1);
        xs = linspace(0,max(H(id)),length(H(id)));
        pv = polyval(pf,xs);
%         plot(xs,pv,line{j},...
%             'color','k',...
%             'linewidth',2),hold on
        yfit = pf(1)*H(id)+pf(2);
        yresid = bdmed(id)-yfit;
        SSresid = sum(yresid.^2);
        SStotal = (length(bdmed(id))-1)*var(bdmed(id));
        rsq(j,i) = 1-SSresid/SStotal;
        fprintf('R-squared: %0.2f\n',rsq(j,i))
    end
end
set(sp,...
    'ylim',[-0.1 0.1],...
    'ytick',-0.1:0.05:0.1,...
    'xlim',[0 1.2],...
    'xtick',0:0.2:1.2)
leg = legend(eb,{'Mudflat';'Fringe';'Forest'});
set(sp(1),'position',[0.13 0.15 0.41 0.76])
set(sp(2),'position',[0.54 0.15 0.41 0.76],...
    'yticklabel',[],'xdir','reverse')
xlabel(sp(1),'Water Depth [m]')
xlabel(sp(2),'Water Depth [m]')
ytext = sprintf('Median elev. relative to\n event beginning [m]');
ylabel(sp(1),ytext)
title(sp(1),'Flood Tide')
title(sp(2),'Ebb Tide')
prettyfigures('text',12,'labels',13,'box',1)
set(leg,'position',[0.87 0.23 0.005 0.005],'box','off')
set(f1,'units','inches');
pos = get(f1,'Position');
% set(f1,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
% print(f1,[sfdir 'Depth_Bdmed'],'-dpdf','-r0')

