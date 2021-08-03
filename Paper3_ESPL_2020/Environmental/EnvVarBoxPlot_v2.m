%Plot boxplots of data fields for manuscript: H, Hs, U^2, eps and taub
%Plot the combined event data
clear,close all
load('d:\Mekong_W2015\DataAnalysis\Paper3\CmbData\CombAllData_flood_v3.mat');
data.flood = dat;
load('d:\Mekong_W2015\DataAnalysis\Paper3\CmbData\CombAllData_ebb_v3.mat');
data.ebb = dat;clear dat
fn = fieldnames(data);
sfdir = 'd:\MemoryStick\GradSchool\DataAnalysis\Paper3\Figures\';
flds = {'Mudflat';'Fringe';'Forest'};
%% All environmental variables (umag, taub, eps)
band = {'flood';'ebb'};
vars = {'umag';'tmax';'eps'};
sp = zeros(3,2);
f1 = figure(1);
set(f1,'PaperOrientation','portrait',...
    'position',[400 100   800   500],...
    'renderer','painters');
w = [1 4;2 5;3 6];
cl = flipud([207 176 126;
    60 166 74;
    4 76 41]./255);
tcr = [0.18,0.26,0.32]; %from Critical_sed_params.m
for i = 1:2
    for j = 1:length(vars)
        sp(w(j,i)) = subplot(2,3,w(j,i));
        x1 = data.(fn{i}).mud.(vars{j})';
        x2 = data.(fn{i}).fringe.(vars{j})';
        x3 = data.(fn{i}).forest.(vars{j})';
        depth = [x1;x2;x3];
        group = [repmat(flds(1),length(x1),1); repmat(flds(2),length(x2),1); repmat(flds(3),length(x3),1)];
        bp = boxplot(depth,group,'symbol','',...
            'colors',[0 0 0],...
            'boxstyle','outline',...
            'datalim',[0 2]);
        set(bp,'linewidth',1.5)
        h = findobj(gca,'Tag','Box');
        for k=1:length(h)
            patch(get(h(k),'XData'),get(h(k),'YData'),cl(k,:),'FaceAlpha',.5);
        end
        if j == 2
            hold on
            group = [1,2,3];
            plot(group,tcr,'*k',...
                'markersize',8,...
                'linewidth',1.5)
            fprintf('Mudflat pct above crit values: %0.2f\n',(nnz(x1>tcr(1))/length(x1))*100)
            fprintf('Fringe pct above crit values: %0.2f\n',(nnz(x2>tcr(2))/length(x2))*100)
            fprintf('Forest pct above crit values: %0.2f\n',(nnz(x3>tcr(3))/length(x3))*100)
        end
    end
end

%Top row
set([sp(1) sp(2) sp(3)],'xticklabel',[])
set(sp(1),'ylim',[0 0.5],'position',[0.07 0.59 0.25 0.35])
set(sp(2),'ylim',[0 2],'position',[0.41 0.59 0.25 0.35])
set(sp(3),'ylim',[0 5E-4],'position',[0.74 0.59 0.25 0.35])
ylabel(sp(1),'$\bar{u} \quad{[m s^{-1}]}$','interpreter','latex')
ylabel(sp(2),'\tau_b [Pa]')
ylabel(sp(3),'\epsilon [Wkg^{-1}]')
%Bottom Row
set(sp(4),'ylim',[0 0.5],'position',[0.07 0.18 0.25 0.35])
set(sp(5),'ylim',[0 2],'position',[0.41 0.18 0.25 0.35])
set(sp(6),'ylim',[0 5E-4],'position',[0.74 0.18 0.25 0.35])
ylabel(sp(4),'$\bar{u} \quad{[m s^{-1}]}$','interpreter','latex')
ylabel(sp(5),'\tau_b [Pa]')
ylabel(sp(6),'\epsilon [Wkg^{-1}]')
prettyfigures('text',12,'labels',13,'box',1)
set(f1,'units','inches');
% pos = get(f1,'Position');
% set(f1,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
% print(f1,[sfdir 'EnvrBoxPlots_v3'],'-dpdf','-r0')
