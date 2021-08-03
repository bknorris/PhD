%Plot bdmed against time, forcing variables against time - GM 2018, Norris et
%al.
clear,close all
load('d:\Projects\Mekong_W2015\DataAnalysis\Paper3\Wavelet\CombEventdata_flood.mat');
data.flood = dat;
load('d:\Projects\Mekong_W2015\DataAnalysis\Paper3\Wavelet\CombEventdata_ebb.mat');
data.ebb = dat;clear dat
fn = fieldnames(data);
sfdir = 'g:\GradSchool\DataAnalysis\Paper3\WorkingFigures\Results\';
flds = {'mud';'fringe';'forest'};

f1 = figure(1);
set(f1,'PaperOrientation','portrait',...
    'position',[500 -50   850   1100],...
    'renderer','opengl');
sp = zeros(4,2);
eb = zeros(3,1);
w = [1 2;3 4;5 6;7 8];
cl = [207 176 126;
    60 166 74;
    4 76 41]./255;
tcr = [0.18,0.26,0.31]; %from Critical_sed_params.m
lspec = {'-';'--';':'};
fields = {'bdmed';'taub';'eps';'umag'};
for i = 1:length(fields)
    for j = 1:length(flds)
        %Flood Tide
        sp(i,1) = subplot(4,2,w(i,1));
        yy = [data.flood.(flds{j}).wave.(fields{i}); data.flood.(flds{j}).ig.(fields{i})];        
        time = [data.flood.(flds{j}).wave.eventl; data.flood.(flds{j}).ig.eventl];
        time = time./80; %max time is roughly 80 mins
%         if j == 3
%             yy(bdmed<-0.06) = [];bdmed(bdmed<-0.06) = [];
%             id = 1:length(yy);
%         else
%             yy(bdmed<-0.1) = [];bdmed(bdmed<-0.1) = [];
%             id = find(bdmed<-0.003 | bdmed > 0.003);
%         end
        bins = linspace(min(time),max(time),15);
        [b,~,q1,q3] = binmedian(time,yy,bins);
        %interp nans
        b(b==0)=NaN;
        q1(q1==0)=NaN;
        if length(find(isnan(q1))) == 15;
            q1 = zeros(15,1);
        end   
        q3(q3==0)=NaN;
        if length(find(isnan(q3))) == 15;
            q3 = zeros(15,1);
        end 
        b = interpnan(b);q1 = interpnan(q1);q3 = interpnan(q3);
        if j == 1
            plot(linspace(0,1,10),zeros(1,10),'-k','linewidth',1.5);hold on
        end
%         if i == 2
%             if j == 2
%                 iid = find(bins > 1E-4,1,'first');
%                 b(iid:end) = b(iid:end)-0.01;q1(iid:end) = q1(iid:end)-0.01;q3(iid:end) = q3(iid:end)-0.01;
%             elseif j == 3
%                 iid = find(bins > 7E-5,1,'first');
%                 b(iid:end) = b(iid:end)-0.01;q1(iid:end) = q1(iid:end)-0.01;q3(iid:end) = q3(iid:end)-0.01;
%             end
%         end
        bins(1) = [];b(1) = [];q1(1) = [];q3(1) = [];
        bl = boundedline(bins,b,[abs(q1) abs(q3)],'alpha','transparency',0.3,...
            'cmap',cl(j,:));hold on
        set(bl,'linewidth',1.5)
        if i == 1,eb(j) = bl;end
        if strcmp(fields{i},'taub')
            for k = 1:3
                plot(0:0.1:1,ones(11,1)*tcr(k),lspec{k},...
                    'linewidth',1,'color','k')
            end
        end
        %Ebb Tide
        sp(i,2) = subplot(4,2,w(i,2));
        yy = [data.ebb.(flds{j}).wave.(fields{i}); data.ebb.(flds{j}).ig.(fields{i})];
        time = [data.ebb.(flds{j}).wave.eventl; data.ebb.(flds{j}).ig.eventl];
        time = time./60; %max time is roughly 60 mins

        bins = linspace(min(time),max(time),15);
        [b,~,q1,q3] = binmedian(time,yy,bins);
        %interp nans
        b(b==0)=NaN;
        q1(q1==0)=NaN;
        if length(find(isnan(q1))) == 15;
            q1 = zeros(15,1);
        end   
        q3(q3==0)=NaN;
        if length(find(isnan(q3))) == 15;
            q3 = zeros(15,1);
        end  
        b = interpnan(b);q1 = interpnan(q1);q3 = interpnan(q3);
        if j == 1
            plot(linspace(0,1,10),zeros(1,10),'-k','linewidth',1.5);hold on
        end
%         if i == 2
%             if j == 2
%                 iid = find(bins > 1E-4,1,'first');
%                 b(iid:end) = b(iid:end)-0.01;q1(iid:end) = q1(iid:end)-0.01;q3(iid:end) = q3(iid:end)-0.01;
%             elseif j == 3
%                 iid = find(bins > 7E-5,1,'first');
%                 b(iid:end) = b(iid:end)-0.01;q1(iid:end) = q1(iid:end)-0.01;q3(iid:end) = q3(iid:end)-0.01;
%             end
%         end
        bins(1) = [];b(1) = [];q1(1) = [];q3(1) = [];
        bl = boundedline(bins,b,[abs(q1) abs(q3)],'alpha','transparency',0.3,...
            'cmap',cl(j,:));hold on
        set(bl,'linewidth',1.5)
        if i == 1,eb(j) = bl;end
        if strcmp(fields{i},'taub')
            for k = 1:3
                plot(0:0.1:1,ones(11,1)*tcr(k),lspec{k},...
                    'linewidth',1,'color','k')
            end
        end
    end
end
%Global Plot adjustments
set(sp,'xlim',[0 1])
set([sp(1,1) sp(1,2)],'ylim',[-0.08 0.04],'ytick',-0.08:0.04:0.04)
set([sp(2,1) sp(2,2)],'ylim',[0 2.5],'ytick',0:0.75:2.5)
set([sp(3,1) sp(3,2)],'ylim',[0 6E-4],'ytick',0:2E-4:6E-4)
set([sp(4,1) sp(4,2)],'ylim',[0 0.4],'ytick',0:0.2:0.4)
set([sp(1,1) sp(1,2) sp(2,1)...
    sp(2,2) sp(3,1) sp(3,2)],...
    'xticklabel',[])
set([sp(1,2) sp(2,2)...
    sp(3,2) sp(4,2)],...
    'yticklabel',[])
set(sp(1,1),'position',[0.12 0.78 0.32 0.18])
set(sp(1,2),'position',[0.5 0.78 0.32 0.18])
set(sp(2,1),'position',[0.12 0.56 0.32 0.18])
set(sp(2,2),'position',[0.5 0.56 0.32 0.18])
set(sp(3,1),'position',[0.12 0.34 0.32 0.18])
set(sp(3,2),'position',[0.5 0.34 0.32 0.18])
set(sp(4,1),'position',[0.12 0.12 0.32 0.18])
set(sp(4,2),'position',[0.5 0.12 0.32 0.18])
leg = legend(eb,{'Mudflat';'Fringe';'Forest'});
txt = sprintf('Median elev. relative\n to initial elev. [m]');
title(sp(1,1),'Flood Tide')
title(sp(1,2),'Ebb Tide')
ylabel(sp(1,1),txt);
ylabel(sp(2,1),'\tau_b [Pa]')
ylabel(sp(3,1),'\epsilon [Wkg^{-1}]')
ylabel(sp(4,1),'$\overline{u} \quad [ms^{-1}]$','interpreter','latex')
xlabel(sp(4,1),'Normalized Event Duration')
xlabel(sp(4,2),'Normalized Event Duration')
prettyfigures('text',12,'labels',13,'box',1)
set(leg,'position',[0.76 0.5 0.005 0.005])
% export_fig([sfdir 'Time_bd_taub_eps_umag'],'-pdf','-nocrop')