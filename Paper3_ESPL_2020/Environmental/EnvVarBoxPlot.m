%Plot boxplots of data fields for manuscript: H, Hs, U^2, eps and taub
%Plot the combined event data
clear,close all
load('e:\Mekong_W2015\DataAnalysis\Paper3\CmbData\CombAllData_flood.mat');
data.flood = dat;
load('e:\Mekong_W2015\DataAnalysis\Paper3\CmbData\CombAllData_ebb.mat');
data.ebb = dat;clear dat
fn = fieldnames(data);
sfdir = 'd:\GradSchool\DataAnalysis\Paper3\Figures\';
flds = {'Mudflat';'Fringe';'Forest'};
%% All environmental variables (H, Hs, Usqd, Eps) excl. Taub
band = {'flood';'ebb'};
vars = {'depth';'Hs';'usqd';'eps'};
sp = zeros(4,2);
f1 = figure(1);
set(f1,'PaperOrientation','portrait',...
    'position',[400 100   1000   500],...
    'renderer','painters');
w = [1 5; 2 6; 3 7;4 8];
cl = flipud([207 176 126;
    60 166 74;
    4 76 41]./255);
for i = 1:2
    for j = 1:length(vars)
        sp(w(j,i)) = subplot(2,4,w(j,i));
        x1 = data.(fn{i}).mud.(vars{j})';
        x2 = data.(fn{i}).fringe.(vars{j})';
        x3 = data.(fn{i}).forest.(vars{j})';
        depth = [x1;x2;x3];
        group = [repmat(flds(1),length(x1),1); repmat(flds(2),length(x2),1); repmat(flds(3),length(x3),1)];
        bp = boxplot(depth,group,'symbol','',...
            'colors',[0 0 0],...
            'boxstyle','outline');
        set(bp,'linewidth',1.5)
        h = findobj(gca,'Tag','Box');
        for k=1:length(h)
            patch(get(h(k),'XData'),get(h(k),'YData'),cl(k,:),'FaceAlpha',.5);
        end
    end
end
%Top row
set(sp(1),'ylim',[0 1.8],'position',[0.07 0.59 0.17 0.35])
set(sp(2),'ylim',[0 0.6],'position',[0.32 0.59 0.17 0.35])
set(sp(3),'ylim',[0 0.04],'position',[0.58 0.59 0.17 0.35])
set(sp(4),'ylim',[0 5E-4],'position',[0.81 0.59 0.17 0.35])
ylabel(sp(1),'Water Depth [m]')
ylabel(sp(2),'Sig. Wave Height [m]')
ylabel(sp(3),'Velocity Squared [m^2s^{-2}]')
ylabel(sp(4),'\epsilon [Wkg^{-1}]')
%Bottom Row
set(sp(5),'ylim',[0 1.8],'position',[0.07 0.12 0.17 0.35])
set(sp(6),'ylim',[0 0.6],'position',[0.32 0.12 0.17 0.35])
set(sp(7),'ylim',[0 0.04],'position',[0.58 0.12 0.17 0.35])
set(sp(8),'ylim',[0 5E-4],'position',[0.81 0.12 0.17 0.35])
ylabel(sp(5),'Water Depth [m]')
ylabel(sp(6),'Sig. Wave Height [m]')
ylabel(sp(7),'Velocity Squared [m^2s^{-2}]')
ylabel(sp(8),'\epsilon [Wkg^{-1}]')
prettyfigures('text',12,'labels',13,'box',1)
set(f1,'units','inches');
pos = get(f1,'Position');
set(f1,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(f1,[sfdir 'EnvrBoxPlots'],'-dpdf','-r0')
%% Taub
f2 = figure(2);
set(f2,'PaperOrientation','portrait',...
    'position',[400 100  300   500],...
    'renderer','painters');
tcr = [0.32,0.4,1.13]; %from Critical_sed_params.m
w = [1 2];
for i = 1:2
    sp(w(1,i)) = subplot(2,1,w(1,i));
    x1 = data.(fn{i}).mud.tmax';
    x2 = data.(fn{i}).fringe.tmax';
    x3 = data.(fn{i}).forest.tmax';
    depth = [x1;x2;x3];
    group = [repmat(flds(1),length(x1),1); repmat(flds(2),length(x2),1); repmat(flds(3),length(x3),1)];
    bp = boxplot(depth,group,'symbol','',...
        'colors',[0 0 0],...
        'boxstyle','outline');
    set(bp,'linewidth',1.5)
    h = findobj(gca,'Tag','Box');
    for k=1:length(h)
        patch(get(h(k),'XData'),get(h(k),'YData'),cl(k,:),'FaceAlpha',.5);
    end
    hold on
    group = [1,2,3];
    plot(group,tcr,'*k',...
        'markersize',8,...
        'linewidth',1.5)
end
set(sp(1),'ylim',[0 2],'position',[0.21 0.58 0.76 0.4])
set(sp(2),'ylim',[0 2],'position',[0.21 0.08 0.76 0.4])
ylabel(sp(1),'\tau_b [Pa]')
ylabel(sp(2),'\tau_b [Pa]')
prettyfigures('text',12,'labels',13,'box',1)
set(f2,'units','inches');
pos = get(f2,'Position');
set(f2,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(f2,[sfdir 'TaubBoxPlots'],'-dpdf','-r0')
