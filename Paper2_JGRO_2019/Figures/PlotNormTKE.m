%Try plotting the u* ratios by canopy height
clear
load('d:\Mekong_W2015\DataAnalysis\Paper2\TKE\NormalizedTKE.mat')
f1 = figure(1);
set(f1,'PaperOrientation','portrait',...
    'position',[400 100 800 600]);
set(gcf,'color','w','paperpositionmode','auto')
symb = {'o','d','p'};symb = repmat(symb,3,1);
sy = repmat({'^'},1,3);symb = [symb;sy];
fn = fieldnames(mydata);
hc = [0.64 0.59 0.61 0.6]; %m, height of canopy
vph = [0.062 0.063 0.061;0.2 0.2 0.2;0.5 0.5 0.5;0.07 0.42 0.81];
zmean = NaN(4,3);
enmean = NaN(4,3);
p = NaN(4,3);
%set colors based on timing of the tide, between low tide and the end of
%the longest VP deployment
t1 = datenum(2015,03,08,13,20,00);
tend = datenum(2015,03,08,18,00,00);
dt = abs(t1-tend);[~,~,~,h,mn,~] = datevec(dt);
telap = h*60+mn;
time = [1 10:10:telap];
toff = [20 40 50 -480];
c = brewermap(length(time),'PuBu');
%%%
for i = 1:4
    dfn = fieldnames(mydata.(fn{i}));
    for ii = 1:3
        vh = vph(i,ii)-0.04-linspace(0.001,0.03,35);
        m = length(mydata.(fn{i}).(dfn{ii}).emax);
        zhc = zeros(m,1);
        for j = 1:length(mydata.(fn{i}).(dfn{ii}).time)
            %color spec
            [~,~,~,h,mn,s] = datevec(mydata.(fn{i}).(dfn{ii}).time(j));
            tstart = h*60+mn+(s/60);
            [~,~,~,h,mn,s] = datevec(t1);
            tidest = h*60+mn+(s/60)+toff(i);
            dt = tstart-tidest;
            tmp = abs(time-dt);
            [~,id] = min(tmp);
            cc = c(id,:);
            %%%
            id = mydata.(fn{i}).(dfn{ii}).eid(j);
            zhc(j) = vh(id)/hc(i);
            p(i,ii) = plot(mydata.(fn{i}).(dfn{ii}).emax(j),zhc(j),...
                symb{i,ii},'color','k',...
                'markersize',8,...
                'markerfacecolor',cc,...
                'linewidth',1);hold on
        end
        zmean(i,ii) = mean(zhc);
        enmean(i,ii) = nanmean(mydata.(fn{i}).(dfn{ii}).emax);
    end
end
%plot mean values
p2(1) = plot(enmean(1:3,1),zmean(1:3,1),'-k',...
    'linewidth',1.5,...
    'marker','.');
p2(2) = plot(enmean(1:3,2),zmean(1:3,2),'--k',...
    'linewidth',1.5,...
    'marker','.');
p2(3) = plot(enmean(1:3,3),zmean(1:3,3),':k',...
    'linewidth',1.5,...
    'marker','.');
p2(4) = plot(enmean(4,:),zmean(4,:),'-.k',...
    'linewidth',1.5,...
    'marker','.');
plot(linspace(0.0001,1,10),ones(1,10)*1,...
    'color','k','linewidth',1)
colormap(c);caxis([0 time(end)])
cb = colorbar('southoutside');
set(cb,'position',[0.13 0.1 0.58 0.03],...
    'xtick',0:50:time(end))
hold off;
set(gca,'ylim',[0 1.4],...
    'ytick',0:0.2:1.4,...
    'xscale','log',...
    'xlim',[10^-4 5*10^-1],...
    'position',...
    [0.13 0.27 0.59 0.7])
tl = xlabel('$\widehat{\epsilon} \quad [-]$');set(tl,'Interpreter','latex','fontname','arial')
cl = xlabel(cb,'Time elapsed after low tide (min)');
ylabel('z/h_c')
leg = legend([p(1,:) p(4,:) p2],{'x = -10 cm';'x = 10 cm';'x = 20 cm';...
    'z/h_c = 0.04';'z/h_c = 0.60';'z/h_c = 1.25';...
    'x = -10 cm';'x = 10 cm';'x = 20 cm';'VTA'});
set(leg,'position',[0.81 0.38 0.1 0.42])
prettyfigures('text',13,'labels',14,'box',1)
savefigdir = 'c:\Users\Bnorr\Documents\GradSchool\DataAnalysis\Paper2\WorkingFigures\VerticalProfiles\';
export_fig([savefigdir 'TKE_prof_HTA_VTA'],'-pdf')
