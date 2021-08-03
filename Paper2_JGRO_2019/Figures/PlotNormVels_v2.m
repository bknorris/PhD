%Try plotting the u,v and w ratios against normalized canopy height
clear
load('d:\Mekong_W2015\DataAnalysis\Paper2\AveragedVelocities\NormalizedVels.mat')
f1 = figure(1);
set(f1,'PaperOrientation','portrait',...
    'position',[400 100   1200   600]);
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
%set colors based on timing of the tide, between low tide and the end of
%the longest VP deployment
t1 = datenum(2015,03,08,13,20,00);
tend = datenum(2015,03,08,19,00,00);
dt = abs(t1-tend);[~,~,~,h,mn,~] = datevec(dt);
telap = h*60+mn;
time = [1 10:10:telap];
toff = [20 40 50 -360];
c = brewermap(length(time),'PuBu');
%%%
for i = 1:4
    dfn = fieldnames(mydata.(fn{i}));        
    for ii = 1:length(dfn)
        vh = vph(i,ii)-0.04-linspace(0.001,0.03,35);
        m = length(mydata.(fn{i}).(dfn{ii}).urav);
        zhc = zeros(m,1);
        for j = 1:m
            %color spec
            [~,~,~,h,mn,s] = datevec(mydata.(fn{i}).(dfn{ii}).utime2(j));
            tstart = h*60+mn+(s/60);
            [~,~,~,h,mn,s] = datevec(t1);
            tidest = h*60+mn+(s/60)+toff(i);
            dt = tstart-tidest;
            tmp = abs(time-dt);
            [~,id] = min(tmp);
            cc = c(id,:);
            %%%
            id = mydata.(fn{i}).(dfn{ii}).uiav(j);
            zhc(j) = vh(id)/hc(i);
            plot(mydata.(fn{i}).(dfn{ii}).urav(j),zhc(j),...
                symb{i,ii},'color','k',...
                'markersize',8,...
                'markerfacecolor',cc,...
                'linewidth',1),hold on
        end
        zmean(i,ii) = mean(zhc);
        urmean(i,ii) = mean(mydata.(fn{i}).(dfn{ii}).urav);
    end
end
%plot mean values
plot(urmean(1:3,1),zmean(1:3,1),'-k',...
    'linewidth',1.5,...
    'marker','.')
plot(urmean(1:3,2),zmean(1:3,2),'--k',...
    'linewidth',1.5,...
    'marker','.')
plot(urmean(1:3,3),zmean(1:3,3),':k',...
    'linewidth',1.5,...
    'marker','.')
plot(urmean(4,:),zmean(4,:),'-.k',...
    'linewidth',1.5,...
    'marker','.')    
hold off;
sp(2) = subplot(132);
zmean = NaN(4,3);
urmean = NaN(4,3);
for i = 1:4
    dfn = fieldnames(mydata.(fn{i}));        
    for ii = 1:length(dfn)
        vh = vph(i,ii)-0.04-linspace(0.001,0.03,35);
        m = length(mydata.(fn{i}).(dfn{ii}).vrav);
        zhc = zeros(m,1);
        for j = 1:m
            %color spec
            [~,~,~,h,mn,s] = datevec(mydata.(fn{i}).(dfn{ii}).vtime2(j));
            tstart = h*60+mn+(s/60);
            [~,~,~,h,mn,s] = datevec(t1);
            tidest = h*60+mn+(s/60)+toff(i);
            dt = tstart-tidest;
            tmp = abs(time-dt);
            [~,id] = min(tmp);
            cc = c(id,:);
            %%%
            id = mydata.(fn{i}).(dfn{ii}).viav(j);
            zhc(j) = vh(id)/hc(i);
            plot(mydata.(fn{i}).(dfn{ii}).vrav(j),zhc(j),...
                symb{i,ii},'color','k',...
                'markersize',8,...
                'markerfacecolor',cc,...
                'linewidth',1),hold on
        end
        zmean(i,ii) = mean(zhc);
        urmean(i,ii) = mean(mydata.(fn{i}).(dfn{ii}).vrav);
    end
end
%plot mean values
plot(urmean(1:3,1),zmean(1:3,1),'-k',...
    'linewidth',1.5,...
    'marker','.')
plot(urmean(1:3,2),zmean(1:3,2),'--k',...
    'linewidth',1.5,...
    'marker','.')
plot(urmean(1:3,3),zmean(1:3,3),':k',...
    'linewidth',1.5,...
    'marker','.')
plot(urmean(4,:),zmean(4,:),'-.k',...
    'linewidth',1.5,...
    'marker','.')    
hold off;
sp(3) = subplot(133);
zmean = NaN(4,3);
urmean = NaN(4,3);
p = NaN(4,3);
for i = 1:4
    dfn = fieldnames(mydata.(fn{i}));        
    for ii = 1:length(dfn)
        vh = vph(i,ii)-0.04-linspace(0.001,0.03,35);
        m = length(mydata.(fn{i}).(dfn{ii}).wrav);
        zhc = zeros(m,1);
        for j = 1:m
            %color spec
            [~,~,~,h,mn,s] = datevec(mydata.(fn{i}).(dfn{ii}).wtime2(j));
            tstart = h*60+mn+(s/60);
            [~,~,~,h,mn,s] = datevec(t1);
            tidest = h*60+mn+(s/60)+toff(i);
            dt = tstart-tidest;
            tmp = abs(time-dt);
            [~,id] = min(tmp);
            cc = c(id,:);
            %%%
            id = mydata.(fn{i}).(dfn{ii}).wiav(j);
            zhc(j) = vh(id)/hc(i);
            p(i,ii) = plot(mydata.(fn{i}).(dfn{ii}).wrav(j),zhc(j),...
                symb{i,ii},'color','k',...
                'markersize',8,...
                'markerfacecolor',cc,...
                'linewidth',1);hold on
        end
        zmean(i,ii) = mean(zhc);
        urmean(i,ii) = mean(mydata.(fn{i}).(dfn{ii}).wrav);
    end
end
%plot mean values
p2(1) = plot(urmean(1:3,1),zmean(1:3,1),'-k',...
    'linewidth',1.5,...
    'marker','.');
p2(2) = plot(urmean(1:3,2),zmean(1:3,2),'--k',...
    'linewidth',1.5,...
    'marker','.');
p2(3) = plot(urmean(1:3,3),zmean(1:3,3),':k',...
    'linewidth',1.5,...
    'marker','.');
p2(4) = plot(urmean(4,:),zmean(4,:),'-.k',...
    'linewidth',1.5,...
    'marker','.');   
colormap(c);caxis([0 time(end)])
cb = colorbar('southoutside');
set(cb,'position',[0.2 0.1 0.5 0.03],...
    'xtick',0:50:time(end))
hold off;
set(sp,'ylim',[0 0.8])
set(sp(1),'xlim',[-2 2],'position',...
    [0.1 0.27 0.21 0.7])
set(sp(2),'xlim',[-2 2],'yticklabel',[],...
    'position',[0.35 0.27 0.21 0.7])
set(sp(3),'xlim',[-10 10],'yticklabel',[],...
    'position',[0.6 0.27 0.21 0.7])
tl = xlabel(sp(1),'$\overline{u}/\overline{u_{h_c}} \quad [-]$');set(tl,'Interpreter','latex','fontname','arial')
tl = xlabel(sp(2),'$\overline{v}/\overline{v_{h_c}} \quad [-]$');set(tl,'Interpreter','latex','fontname','arial')
tl = xlabel(sp(3),'$\overline{w}/\overline{w_{h_c}} \quad [-]$');set(tl,'Interpreter','latex','fontname','arial')
cl = xlabel(cb,'Time elapsed after low tide (min)');
ylabel(sp(1),'z/h_c')
leg = legend([p(1,:) p(4,1:2) p2],{'x = -10 cm';'x = 10 cm';'x = 20 cm';'z/h_c = 0.04';'z/h_c = 0.60';...
    'x = -10 cm';'x = 10 cm';'x = 20 cm';'VTA'});
set(leg,'position',[0.85 0.38 0.1 0.42])
prettyfigures('text',13,'labels',14,'box',1)
savefigdir = 'c:\Users\Bnorr\Documents\GradSchool\DataAnalysis\Paper2\WorkingFigures\VerticalProfiles\';
export_fig([savefigdir 'Vel_prof_HTA_VTA'],'-pdf')