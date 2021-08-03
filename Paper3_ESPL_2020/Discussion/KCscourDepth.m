%Calculate KC numbers, compare to max(erosion)
clear,close all
load('d:\Mekong_W2015\DataAnalysis\Paper3\CmbData\CombAllData_flood.mat');
data.flood = dat;
load('d:\Mekong_W2015\DataAnalysis\Paper3\CmbData\CombAllData_ebb.mat');
data.ebb = dat;clear dat
fn = fieldnames(data);
sfdir = 'd:\GradSchool\DataAnalysis\Paper3\Figures\';
flds = {'fringe';'forest'};
for i = 1:2
    kc = zeros(2,6);
    sd = zeros(2,6);
    for j = 1:length(flds)
        exps = unique(data.(fn{i}).(flds{j}).expnum,'stable');
        d = unique(data.(fn{i}).(flds{j}).d,'stable').*1.5;
        dbd = data.(fn{i}).(flds{j}).bdist;
        KC = data.(fn{i}).(flds{j}).KC(1:length(dbd));
        expn = data.(fn{i}).(flds{j}).expnum(1:length(dbd));
        if i == i && j == 1
            exps(5) = [];
        elseif i == 2 && j == 1
            exps(3) = [];
        end
        for k = 1:length(exps)
            idx = find(expn == exps(k));
            bd = dbd(idx);
            kc(j,k) = max(KC(idx));
            sd(j,k) = abs(min(bd))/d(k);
        end
    end
    data.(fn{i}).kc = kc;
    data.(fn{i}).sd = sd;
end
f1 = figure(1);
set(f1,'PaperOrientation','portrait',...
    'position',[400 100   500   400],...
    'renderer','painters');hold on
symb = {'^';'s'};
line = {'-.';'--'};
cl = [0.7 0.7 0.7;0.3 0.3 0.3];
pi = zeros(2,1);
for i = 1:2
    kc = data.(fn{i}).kc;
    xs = [kc(1,:) kc(2,:)];xs(xs == 0) = [];
    sd = data.(fn{i}).sd;
    ys = [sd(1,:) sd(2,:)];ys(ys == 0) = [];
    %calculate exponential polyfit
    pf = polyfit(xs,ys,2);
    xx = linspace(1,150,1000);
    if i == 1
        sfkc = 1.3*(1-exp(-0.03*(xx-6)));
        sf = plot(xx,sfkc,'-',...
            'color','k',...
            'linewidth',1.5);
    end
    for j = 1:2
        a = kc(j,:);a(a == 0) = [];
        b = sd(j,:);b(b == 0) = [];
        pi(j) = plot(a,b,symb{j},...
            'markerfacecolor',cl(i,:),...
            'markeredgecolor','k',...
            'markersize',6);
    end
end
leg = legend([pi; sf],{'Fringe';'Forest';'SF1992'});
set(leg,'position',[0.7 0.2 0.05 0.05])
xlabel('KC')
ylabel('S/D')
set(gca,...
    'xscale','log',...
    'yscale','log',...
    'xlim',[10 150],...
    'ylim',[0.1 5])
prettyfigures('text',12,'labels',13,'box',1)
set(f1,'units','inches');
pos = get(f1,'Position');
set(f1,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
% print(f1,[sfdir 'KC_sdepth'],'-dpdf','-r0')
