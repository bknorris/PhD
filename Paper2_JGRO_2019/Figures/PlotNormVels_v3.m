%Try plotting the u,v and w ratios against normalized canopy height
clear, close all
load('d:\Mekong_W2015\DataAnalysis\Paper2\AveragedVelocities\NormalizedVels.mat')
f1 = figure(1);
set(f1,'PaperOrientation','portrait',...
    'position',[400 100   1200   400]);
set(gcf,'color','w','paperpositionmode','auto')
symb = {'o','d','p'};symb = repmat(symb,3,1);
sy = repmat({'^'},1,3);symb = [symb;sy];
c = flipud([0.2 0.2 0.2;0.5 0.5 0.5;0.7 0.7 0.7]);
fn = fieldnames(mydata);
hold on
hc = [0.64 0.59 0.61 0.6]; %m, height of canopy
vph = [0.062 0.063 0.061;0.2 0.2 0.2;0.5 0.5 0.5;0.07 0.42 0];
sp(1) = subplot(131);
zmean = NaN(4,3);
urmean = NaN(4,3);
cmap = csvread('c:\Users\Bnorr\Documents\GradSchool\Design\Qualitative8_4.csv');
cmap = flipud(cmap)./255;
for i = 1:4
    dfn = fieldnames(mydata.(fn{i}));
    for ii = 1:length(dfn)
        if i < 4
            cl = cmap(ii,:);
        else
            cl = cmap(i,:);
        end
        vh = vph(i,ii)-0.04-linspace(0.001,0.03,35);
        id = round(mean(mydata.(fn{i}).(dfn{ii}).uiav));
        xs = nanmean(mydata.(fn{i}).(dfn{ii}).urav);
        xstd = std(mydata.(fn{i}).(dfn{ii}).urav,1);
        zhc = vh(id)/hc(i);
        eb = ploterr(xs,zhc,xstd,[],'abshhx',0.05);set(eb,'color',cl,'linewidth',1.5),hold on
        zmean(i,ii) = mean(zhc);
        urmean(i,ii) = xs;
    end  
end
%plot mean values
plot(urmean(1:3,1),zmean(1:3,1),'-',...
    'color',cmap(1,:),...
    'linewidth',1.5,...
    'marker','.')
plot(urmean(1:3,2),zmean(1:3,2),'--',...
    'color',cmap(2,:),...
    'linewidth',1.5,...
    'marker','.')
plot(urmean(1:3,3),zmean(1:3,3),':',...
    'color',cmap(3,:),...
    'linewidth',1.5,...
    'marker','.')
plot(urmean(4,:),zmean(4,:),'-.',...
    'color',cmap(4,:),...
    'linewidth',1.5,...
    'marker','.')
for i = 1:4
    dfn = fieldnames(mydata.(fn{i}));
    for ii = 1:length(dfn)
        if i < 4
            cl = cmap(ii,:);
        else
            cl = cmap(i,:);
        end
        vh = vph(i,ii)-0.04-linspace(0.001,0.03,35);
        id = round(mean(mydata.(fn{i}).(dfn{ii}).uiav));
        xs = nanmean(mydata.(fn{i}).(dfn{ii}).urav);
        zhc = vh(id)/hc(i);
        plot(xs,zhc,...
            symb{i,ii},'color','k',...
            'markersize',9,...
            'markerfacecolor',cl,...
            'linewidth',1),hold on
    end
end
hold off;
sp(2) = subplot(132);
zmean = NaN(4,3);
vrmean = NaN(4,3);
for i = 1:4
    dfn = fieldnames(mydata.(fn{i}));
    for ii = 1:length(dfn)
        if i < 4
            cl = cmap(ii,:);
        else
            cl = cmap(i,:);
        end
        vh = vph(i,ii)-0.04-linspace(0.001,0.03,35);
        id = round(mean(mydata.(fn{i}).(dfn{ii}).viav));
        xs = nanmean(mydata.(fn{i}).(dfn{ii}).vrav);
        xstd = std(mydata.(fn{i}).(dfn{ii}).vrav,1);
        zhc = vh(id)/hc(i);
        eb = ploterr(xs,zhc,xstd,[],'abshhx',0.05);set(eb,'color',cl,'linewidth',1.5),hold on
        zmean(i,ii) = mean(zhc);
        vrmean(i,ii) = xs;
    end  
end
%plot mean values
plot(vrmean(1:3,1),zmean(1:3,1),'-',...
    'color',cmap(1,:),...
    'linewidth',1.5,...
    'marker','.')
plot(vrmean(1:3,2),zmean(1:3,2),'--',...
    'color',cmap(2,:),...
    'linewidth',1.5,...
    'marker','.')
plot(vrmean(1:3,3),zmean(1:3,3),':',...
    'color',cmap(3,:),...
    'linewidth',1.5,...
    'marker','.')
plot(vrmean(4,:),zmean(4,:),'-.',...
    'color',cmap(4,:),...
    'linewidth',1.5,...
    'marker','.')
for i = 1:4
    dfn = fieldnames(mydata.(fn{i}));
    for ii = 1:length(dfn)
        if i < 4
            cl = cmap(ii,:);
        else
            cl = cmap(i,:);
        end
        vh = vph(i,ii)-0.04-linspace(0.001,0.03,35);
        id = round(mean(mydata.(fn{i}).(dfn{ii}).viav));
        xs = nanmean(mydata.(fn{i}).(dfn{ii}).vrav);
        zhc = vh(id)/hc(i);
        plot(xs,zhc,...
            symb{i,ii},'color','k',...
            'markersize',9,...
            'markerfacecolor',cl,...
            'linewidth',1),hold on
    end
end
hold off;
sp(3) = subplot(133);
zmean = NaN(4,3);
wrmean = NaN(4,3);
for i = 1:4
    dfn = fieldnames(mydata.(fn{i}));
    for ii = 1:length(dfn)
        if i < 4
            cl = cmap(ii,:);
        else
            cl = cmap(i,:);
        end
        vh = vph(i,ii)-0.04-linspace(0.001,0.03,35);
        id = round(mean(mydata.(fn{i}).(dfn{ii}).wiav));
        xs = nanmean(mydata.(fn{i}).(dfn{ii}).wrav);
        xstd = std(mydata.(fn{i}).(dfn{ii}).wrav,1);
        zhc = vh(id)/hc(i);
        eb = ploterr(xs,zhc,xstd,[],'abshhx',0.05);set(eb,'color',cl,'linewidth',1.5),hold on
        zmean(i,ii) = mean(zhc);
        wrmean(i,ii) = xs;
    end  
end
%plot mean values
p2(1) = plot(wrmean(1:3,1),zmean(1:3,1),'-',...
    'color',cmap(1,:),...
    'linewidth',1.5,...
    'marker','.');
p2(2) = plot(wrmean(1:3,2),zmean(1:3,2),'--',...
    'color',cmap(2,:),...
    'linewidth',1.5,...
    'marker','.');
p2(3) = plot(wrmean(1:3,3),zmean(1:3,3),':',...
    'color',cmap(3,:),...
    'linewidth',1.5,...
    'marker','.');
p2(4) = plot(wrmean(4,:),zmean(4,:),'-.',...
    'color',cmap(4,:),...
    'linewidth',1.5,...
    'marker','.');
p = NaN(4,3);
for i = 1:4
    dfn = fieldnames(mydata.(fn{i}));
    for ii = 1:length(dfn)
        if i < 4
            cl = cmap(ii,:);
        else
            cl = cmap(i,:);
        end
        vh = vph(i,ii)-0.04-linspace(0.001,0.03,35);
        id = round(mean(mydata.(fn{i}).(dfn{ii}).wiav));
        xs = nanmean(mydata.(fn{i}).(dfn{ii}).wrav);
        zhc = vh(id)/hc(i);
        p(i,ii) = plot(xs,zhc,...
            symb{i,ii},'color','k',...
            'markersize',9,...
            'markerfacecolor',cl,...
            'linewidth',1);hold on
    end
end
hold off;
set(sp,'ylim',[0 0.8])
set(sp(1),'xlim',[-2 2],...
    'xtick',-2:1:2,...
    'position',...
    [0.1 0.17 0.21 0.75])
set(sp(2),'xlim',[-2 2],...
    'xtick',-2:1:2,...
    'yticklabel',[],...
    'position',[0.35 0.17 0.21 0.75])
set(sp(3),'xlim',[-10 10],...
    'xtick',-10:5:10,...
    'yticklabel',[],...
    'position',[0.6 0.17 0.21 0.75])
tl = xlabel(sp(1),'$\overline{u}/\overline{u_{h_c}} \quad [-]$');set(tl,'Interpreter','latex','fontname','arial')
tl = xlabel(sp(2),'$\overline{v}/\overline{v_{h_c}} \quad [-]$');set(tl,'Interpreter','latex','fontname','arial')
tl = xlabel(sp(3),'$\overline{w}/\overline{w_{h_c}} \quad [-]$');set(tl,'Interpreter','latex','fontname','arial')
ylabel(sp(1),'z/h_c')
leg = legend([p(1,:) p(4,1:2) p2],{'x = -10 cm';'x = 10 cm';'x = 20 cm';'z/h_c = 0.04';'z/h_c = 0.60';...
    'x = -10 cm';'x = 10 cm';'x = 20 cm';'VTA'});
set(leg,'position',[0.85 0.34 0.1 0.42])
prettyfigures('text',13,'labels',14,'box',1,'grid',1,'gstyle',':')
savefigdir = 'c:\Users\Bnorr\Documents\GradSchool\DataAnalysis\Paper2\WorkingFigures\VerticalProfiles\';
export_fig([savefigdir 'Vel_prof_HTA_VTA_v2'],'-pdf')