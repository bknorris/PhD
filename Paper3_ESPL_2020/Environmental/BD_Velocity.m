%Plot bdmed against forcing variables - GM 2018, Norris et
%al.
clear,close all
load('d:\Mekong_W2015\DataAnalysis\Paper3\Wavelet\CombEventdata_flood.mat');
data.flood = dat;
load('d:\Mekong_W2015\DataAnalysis\Paper3\Wavelet\CombEventdata_ebb.mat');
data.ebb = dat;clear dat
fn = fieldnames(data);
sfdir = 'e:\GradSchool\DataAnalysis\Paper3\WorkingFigures\Results\';
flds = {'mud';'fringe';'forest'};

f1 = figure(1);
set(f1,'PaperOrientation','portrait',...
    'position',[400 100   400   500],...
    'renderer','painters');
w = [1 2]; %window order
sp = zeros(1,2);
eb = zeros(3,1);
cl = [207 176 126;
    60 166 74;
    4 76 41]./255;
Ucr = [0.18,0.22,0.20]; %from Critical_sed_params.m
lspec = {'-';'--';':'};
fields = {'umag'};
for i = 1:length(fields)
    for j = 1:length(flds)
        %Flood Tide
        sp(w(i,1)) = subplot(2,1,w(i,1));
        bdmed = [data.flood.(flds{j}).wave.bdmed; data.flood.(flds{j}).ig.bdmed];
        xx = [data.flood.(flds{j}).wave.(fields{i}); data.flood.(flds{j}).ig.(fields{i})];
        if j == 3
            xx(bdmed<-0.06) = [];bdmed(bdmed<-0.06) = [];
            id = 1:length(xx);
        else
            xx(bdmed<-0.1) = [];bdmed(bdmed<-0.1) = [];
            id = find(bdmed<-0.003 | bdmed > 0.003);
        end
        bins = linspace(min(xx(id)),max(xx(id)),15);
        [b,~,q1,q3] = binmedian(xx(id),bdmed(id),bins);
        %interp nans
        q1(q1==0)=NaN;q3(q3==0)=NaN;
        b = interpnan(b);q1 = interpnan(q1);q3 = interpnan(q3);
        if j == 1
            plot(linspace(0,2,10),zeros(1,10),'--k','linewidth',1.5);hold on
        end
        bins(1) = [];b(1) = [];q1(1) = [];q3(1) = [];
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
        sp(w(i,2)) = subplot(2,1,w(i,2));
        bdmed = [data.ebb.(flds{j}).wave.bdmed; data.ebb.(flds{j}).ig.bdmed];
        xx = [data.ebb.(flds{j}).wave.(fields{i}); data.ebb.(flds{j}).ig.(fields{i})];
        if j == 1
            plot(linspace(0,2,10),zeros(1,10),'--k','linewidth',1.5);hold on
        end
        if j == 3
            xx(bdmed<-0.06) = [];bdmed(bdmed<-0.06) = [];
            id = 1:length(xx);
        else
            xx(bdmed<-0.1) = [];bdmed(bdmed<-0.1) = [];
            id = find(bdmed<-0.003 | bdmed > 0.003);
        end
        bins = linspace(min(xx(id)),max(xx(id)),15);
        [b,~,q1,q3] = binmedian(xx(id),bdmed(id),bins);
        %interp nans
        q1(q1==0)=NaN;q3(q3==0)=NaN;
        b = interpnan(b);q1 = interpnan(q1);q3 = interpnan(q3);
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
% Global Plot adjustments
set(sp,...
    'xlim',[0 0.52],...
    'xtick',0:0.1:0.5,...
    'ylim',[-0.1 0.1])
set(sp(1),...
    'xticklabel',[])
%Positioning
set(sp(1),'position',[0.18 0.57 0.7 0.4])
set(sp(2),'position',[0.18 0.12 0.7 0.4])
%Labeling
leg = legend(eb,{'Mudflat';'Fringe';'Forest'});
suplabel('Median elev. relative to event beginning [m]','y');
xlabel(sp(2),'Mean Current Velocity [ms^{-1}]')
prettyfigures('text',12,'labels',13,'box',1)
set(leg,'position',[0.75 0.2 0.005 0.005],'box','off')
set(f1,'units','inches');
pos = get(f1,'Position');
set(f1,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(f1,[sfdir 'Bdmed_Umag'],'-dpdf','-r0')
