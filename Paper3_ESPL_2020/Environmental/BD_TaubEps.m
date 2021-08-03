%Plot bdmed against forcing variables - GM 2018, Norris et
%al.
clear,close all
load('e:\Mekong_W2015\DataAnalysis\Paper3\Wavelet\CombEventdata_flood.mat');
data.flood = dat;
load('e:\Mekong_W2015\DataAnalysis\Paper3\Wavelet\CombEventdata_ebb.mat');
data.ebb = dat;clear dat
fn = fieldnames(data);
sfdir = 'f:\GradSchool\DataAnalysis\Paper3\WorkingFigures\Results\';
flds = {'mud';'fringe';'forest'};

f1 = figure(1);
set(f1,'PaperOrientation','portrait',...
    'position',[400 100   800   500],...
    'renderer','painters');
w = [1 3;2 4]; %window order
sp = zeros(2,2);
eb = zeros(3,1);
cl = [207 176 126;
    60 166 74;
    4 76 41]./255;
tcr = [0.18,0.26,0.31]; %from Critical_sed_params.m
lspec = {'-';'--';':'};
fields = {'taub';'eps'};
for i = 1:length(fields)
    for j = 1:length(flds)
        %Flood Tide
        sp(w(i,1)) = subplot(2,2,w(i,1));
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
            plot(linspace(1E-6,100,10),zeros(1,10),'-k','linewidth',1.5);hold on
        end
        if i == 2
            if j == 2
                iid = find(bins > 1E-4,1,'first');
                b(iid:end) = b(iid:end)-0.01;q1(iid:end) = q1(iid:end)-0.01;q3(iid:end) = q3(iid:end)-0.01;
            elseif j == 3
                iid = find(bins > 7E-5,1,'first');
                b(iid:end) = b(iid:end)-0.01;q1(iid:end) = q1(iid:end)-0.01;q3(iid:end) = q3(iid:end)-0.01;
            end
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
        sp(w(i,2)) = subplot(2,2,w(i,2));
        bdmed = [data.ebb.(flds{j}).wave.bdmed; data.ebb.(flds{j}).ig.bdmed];
        xx = [data.ebb.(flds{j}).wave.(fields{i}); data.ebb.(flds{j}).ig.(fields{i})];
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
            plot(linspace(1E-6,100,10),zeros(1,10),'-k','linewidth',1.5);hold on
        end
        if i == 1
            if j == 3
                iid = find(bins < 0.13,1,'last');
                b(1:iid) = b(1:iid)+0.002;q1(1:iid) = q1(1:iid)+0.002;q3(1:iid) = q3(1:iid)+0.002;
            end
        end
        if i == 2
            if j == 3
                iid = find(bins < 7E-5,1,'last');
                b(1:iid) = b(1:iid)+0.01;q1(1:iid) = q1(1:iid)+0.01;q3(1:iid) = q3(1:iid)+0.01;
            end
        end
        bins(1) = [];b(1) = [];q1(1) = [];q3(1) = [];
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
%Global Plot adjustments
set([sp(1) sp(3)],...
    'xlim',[0 1.5],...
    'ylim',[-0.1 0.1])
set([sp(2) sp(4)],...
    'xlim',[8E-6 5E-4],...
    'ylim',[-0.1 0.1],...
    'xscale','log')
set(sp(1),...
    'xticklabel',[])
set(sp(2),...
    'yticklabel',[],...
    'xticklabel',[])
set(sp(4),...
    'yticklabel',[])
%Positioning
set(sp(1),'position',[0.12 0.57 0.38 0.4])
set(sp(2),'position',[0.54 0.57 0.38 0.4])
set(sp(3),'position',[0.12 0.12 0.38 0.4])
set(sp(4),'position',[0.54 0.12 0.38 0.4])
%Labeling
leg = legend(eb,{'Mudflat';'Fringe';'Forest'});
suplabel('Median elev. relative to event beginning [m]','y');
xlabel(sp(3),'\tau_b [Pa]')
xlabel(sp(4),'\epsilon [Wkg^{-1}]')
prettyfigures('text',12,'labels',13,'box',1)
set(leg,'position',[0.85 0.63 0.005 0.005],'box','off')
set(f1,'units','inches');
pos = get(f1,'Position');
% set(f1,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
% print(f1,[sfdir 'Bdmed_TaubEps'],'-dpdf','-r0')
