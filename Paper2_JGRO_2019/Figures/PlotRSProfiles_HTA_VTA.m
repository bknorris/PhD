%Plot averaged Reynolds stresses for the HTA and VTA (2x3 panel subplots).
%This is version 2 of this script.

clear
load('D:\Projects\Mekong_W2015\DataAnalysis\Paper2\RS_estimates\RStress_10min_v4test.mat')
savefigdir = 'f:\GradSchool\DataAnalysis\Paper2\WorkingFigures\VerticalProfiles\';

%plot routine
f1 = figure;
set(f1,'PaperOrientation','portrait',...
    'position',[400 100   600  800]);
set(gcf,'color','w','paperpositionmode','auto')
s = [5 3 1]; %HTA
p = [6 4 2]; %VTA
sp1 = zeros(3,1); %HTA
sp2 = zeros(3,1); %VTA
dn = fieldnames(RS);
hab = [0.066 0.066 0.066;
    0.248 0.248 0.248;
    0.550 0.550 0.550;
    0.06 0.416 0.806];
hc = [0.58 0.58 0.58 0.6];
symb = {'o';'d';'p'};
line = {'-';'--';'-.'};
cl = flipud([0.1 0.1 0.1;0.5 0.5 0.5;0.7 0.7 0.7]);
for i = 1:4
    fn = fieldnames(RS.(dn{i}));
    for j = 1:3
        tmp = RS.(dn{i}).(fn{j}).uw;
        if i == 1
            tmp = tmp(25:80,:);
        end
        if i == 4
            tmp = tmp(15:40,:);
        end
        id = find(tmp> 0.25 | tmp < -0.25);
        tmp(id) = NaN;
%         if i == 1 && j == 2
%             uw = nanmean(tmp)./-30;
%         else
%             uw = nanmean(tmp)./10;
%         end
        uw = nanmean(tmp);
        z = hab(i,j)-0.04-linspace(0,0.03,35);
        zhc = z./hc(i);
        markx = uw(1:5:end);
        marky = zhc(1:5:end);
        if i < 4
            sp1(i) = subplot(3,2,s(i));
            plot(uw,zhc,line{j},...
                'color',cl(j,:),...
                'linewidth',1.5), hold on
            plot(markx,marky,symb{j},...
                'linewidth',1.5,...
                'markerfacecolor',cl(j,:),...
                'markeredgecolor','k',...
                'markersize',6), hold on
            grid on
        else
            sp2(j) = subplot(3,2,p(j));
            plot(uw,zhc,'-',...
                'color','k',...
                'linewidth',1.5), hold on
            plot(markx,marky,'^',...
                'linewidth',1.5,...
                'markerfacecolor','w',...
                'markeredgecolor','k',...
                'markersize',6)
            grid on
        end
    end
end
break
%global plot adjustments
set(sp1,'xlim',[-0.02 0.03],...
    'xtick',-0.02:0.02:0.04,'gridlinestyle',':')
set(sp2,'xlim',[-0.005 0.005],...
    'xtick',-0.005:0.0025:0.005,'gridlinestyle',':')
set([sp1(2) sp1(3)],'xticklabel',[])
set([sp2(2) sp2(3)],'xticklabel',[])
%HTA
set(sp1(1),'ylim',[0 0.05],...
    'ytick',0:0.015:0.05,...
    'position',[0.13 0.1 0.35 0.25])
set(sp1(2),'ylim',[0.3 0.36],...
    'ytick',0.29:0.02:0.35,...
    'position',[0.13 0.4 0.35 0.25])
set(sp1(3),'ylim',[0.82 0.88],...
    'ytick',0.82:0.02:0.88,...
    'position',[0.13 0.7 0.35 0.25])
title(sp1(3),'HTA')
ylabel(sp1(2),'z/h_c')
xl = xlabel(sp1(1),'$-\overline{u''w''} \quad (m^2s^{-2})$');
set(xl,'Interpreter','latex','fontname','arial')
%VTA
set(sp2(1),'ylim',[0.01 0.05],...
    'ytick',0.01:0.01:0.05,...
    'position',[0.6 0.1 0.35 0.25])
set(sp2(2),'ylim',[0.59 0.63],...
    'ytick',0.58:0.015:0.63,...
    'position',[0.6 0.4 0.35 0.25])
set(sp2(3),'ylim',[1.22 1.28],...
    'ytick',1.22:0.02:1.28,...
    'position',[0.6 0.7 0.35 0.25])
title(sp2(3),'VTA')
xl = xlabel(sp2(1),'$-\overline{u''w''} \quad (m^2s^{-2})$');
set(xl,'Interpreter','latex','fontname','arial')
prettyfigures('text',13,'labels',14,'box',1,'tlength',[0.025 0.025])
% export_fig([savefigdir 'RSprof_HTA_VTA'],'-pdf')
