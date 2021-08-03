%Plot normalized eddy viscosity
clear, close all
load('d:\Mekong_W2015\DataAnalysis\Paper2\RS_estimates\NormEddyVis.mat')
f1 = figure(1);
set(f1,'PaperOrientation','portrait',...
    'position',[400 100 700 500]);
set(gcf,'color','w','paperpositionmode','auto')
symb = {'o','d','p'};symb = repmat(symb,3,1);
sy = repmat({'^'},1,3);symb = [symb;sy];
fn = fieldnames(mydata);
hc = [0.64 0.59 0.61 0.6]; %m, height of canopy
vph = [0.062 0.063 0.061;0.2 0.2 0.2;0.5 0.5 0.5;0.07 0.42 0.81];
zmean = NaN(4,3);
nmean = NaN(4,3);
p = NaN(4,3);
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
        id = round(mean(mydata.(fn{i}).(dfn{ii}).nid));
        xs = nanmean(mydata.(fn{i}).(dfn{ii}).nmax);
        xstd = nanstd(mydata.(fn{i}).(dfn{ii}).nmax);
        if xs-xstd < 0
            xsd = xs-abs(xstd-xs);
        else
            xsd = [xstd xstd];
        end
        zhc = vh(id)/hc(i);
        eb = ploterr(xs,zhc,xsd,[],'abshhx',0.05);set(eb,'color',cl,'linewidth',1.5),hold on
        zmean(i,ii) = mean(zhc);
        nmean(i,ii) = nanmean(mydata.(fn{i}).(dfn{ii}).nmax);
    end
end
%plot mean values
p2(1) = plot(nmean(1:3,1),zmean(1:3,1),'-',...
    'color',cmap(1,:),...
    'linewidth',1.5,...
    'marker','.');
p2(2) = plot(nmean(1:3,2),zmean(1:3,2),'--',...
    'color',cmap(2,:),...
    'linewidth',1.5,...
    'marker','.');
p2(3) = plot(nmean(1:3,3),zmean(1:3,3),':',...
    'color',cmap(3,:),...
    'linewidth',1.5,...
    'marker','.');
p2(4) = plot(nmean(4,:),zmean(4,:),'-.',...
    'color',cmap(4,:),...
    'linewidth',1.5,...
    'marker','.');
for i = 1:4
    dfn = fieldnames(mydata.(fn{i}));
    for ii = 1:length(dfn)
        if i < 4
            cl = cmap(ii,:);
        else
            cl = cmap(i,:);
        end
        vh = vph(i,ii)-0.04-linspace(0.001,0.03,35);
        id = round(mean(mydata.(fn{i}).(dfn{ii}).nid));
        xs = nanmean(mydata.(fn{i}).(dfn{ii}).nmax);
        zhc = vh(id)/hc(i);
        p(i,ii) = plot(xs,zhc,...
            symb{i,ii},'color','k',...
            'markersize',8,...
            'markerfacecolor',cl,...
            'linewidth',1);
    end
end

