%Try plotting the u* ratios by canopy height
clear
load('d:\Mekong_W2015\DataAnalysis\Paper2\TKE\NormalizedTKE.mat')
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
enmean = NaN(4,3);
cmap = csvread('c:\Users\Bnorr\Documents\GradSchool\Design\Qualitative8_4.csv');
cmap = flipud(cmap)./255;
for i = 1:4
    dfn = fieldnames(mydata.(fn{i}));
    for ii = 1:3
        if i < 4
            cl = cmap(ii,:);
        else
            cl = cmap(i,:);
        end
        vh = vph(i,ii)-0.04-linspace(0.001,0.03,35);
        id = round(mean(mydata.(fn{i}).(dfn{ii}).eid));
        zhc = vh(id)/hc(i);
        xs = nanmean(mydata.(fn{i}).(dfn{ii}).emax);
        xstd = nanstd(mydata.(fn{i}).(dfn{ii}).emax,1);
        if xs-xstd < 0
            xsd = xs-abs(xstd-xs);
        else
            xsd = [xstd xstd];
        end
        eb = ploterr(xs,zhc,xsd,[],'abshhx',0.05);set(eb,'color',cl,'linewidth',1.5),hold on
        zmean(i,ii) = mean(zhc);
        enmean(i,ii) = nanmean(mydata.(fn{i}).(dfn{ii}).emax);
    end
end
%plot mean values
p2(1) = plot(enmean(1:3,1),zmean(1:3,1),'-',...
    'color',cmap(1,:),...
    'linewidth',1.5,...
    'marker','.');
p2(2) = plot(enmean(1:3,2),zmean(1:3,2),'--',...
    'color',cmap(2,:),...
    'linewidth',1.5,...
    'marker','.');
p2(3) = plot(enmean(1:3,3),zmean(1:3,3),':',...
    'color',cmap(3,:),...
    'linewidth',1.5,...
    'marker','.');
p2(4) = plot(enmean(4,:),zmean(4,:),'-.',...
    'color',cmap(4,:),...
    'linewidth',1.5,...
    'marker','.');
plot(linspace(0.0001,1,10),ones(1,10)*1,...
    'color','k','linewidth',1)
p = NaN(4,3);
for i = 1:4
    dfn = fieldnames(mydata.(fn{i}));
    for ii = 1:3
        if i < 4
            cl = cmap(ii,:);
        else
            cl = cmap(i,:);
        end
        vh = vph(i,ii)-0.04-linspace(0.001,0.03,35);
        id = round(mean(mydata.(fn{i}).(dfn{ii}).eid));
        zhc = vh(id)/hc(i);
        xs = nanmean(mydata.(fn{i}).(dfn{ii}).emax);
        p(i,ii) = plot(xs,zhc,...
            symb{i,ii},'color','k',...
            'markersize',8,...
            'markerfacecolor',cl,...
            'linewidth',1);hold on
    end
end
hold off;
set(gca,'ylim',[0 1.4],...
    'ytick',0:0.2:1.4,...
    'xscale','log',...
    'position',...
    [0.13 0.17 0.58 0.76])
tl = xlabel('$\widehat{\epsilon} \quad [-]$');set(tl,'Interpreter','latex','fontname','arial')
ylabel('z/h_c')
leg = legend([p(1,:) p(4,:) p2],{'x = -10 cm';'x = 10 cm';'x = 20 cm';...
    'z/h_c = 0.04';'z/h_c = 0.60';'z/h_c = 1.25';...
    'x = -10 cm';'x = 10 cm';'x = 20 cm';'VTA'});
set(leg,'position',[0.81 0.34 0.1 0.42])
prettyfigures('text',13,'labels',14,'box',1)
savefigdir = 'c:\Users\Bnorr\Documents\GradSchool\DataAnalysis\Paper2\WorkingFigures\VerticalProfiles\';
export_fig([savefigdir 'TKE_prof_HTA_VTA_v2'],'-pdf')
