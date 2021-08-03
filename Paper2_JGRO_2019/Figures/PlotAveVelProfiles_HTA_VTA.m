%Plot time-averaged velocity profiles for the HTA and VTA experiments. This
%is version 2 of this script.

clear
datdir = 'd:\Projects\Mekong_W2015\DataAnalysis\Paper2\AveragedVelocities\';
files = {'HTA_1Vave.mat';'HTA_2Vave.mat';'HTA_4Vave.mat';...
    'VTA_2vp1Vave.mat';'VTA_2vp2Vave.mat';'VTA_2vp3Vave.mat'};
savefigdir = 'd:\Projects\Mekong_W2015\Figures\Paper1\';

%plot routine
f1 = figure;
set(f1,'PaperOrientation','portrait',...
    'position',[400 100   600  800]);
set(gcf,'color','w','paperpositionmode','auto')
s = [5 3 1 6 4 2]; %subplot order
sp = zeros(6,1);
hab = [0.062 0.063 0.061;
    0.242 0.240 0.240;
    0.550 0.550 0.550;
    0.07  0 0; %zeros are dummy to get the indexing right with 6 files
    0.416 0 0;
    0.806 0 0];
hc = [0.58 0.58 0.58 0.6 0.6 0.6];
symb = {'o';'d';'p'};
line = {'-';'--';'-.'};
cl = flipud([0.1 0.1 0.1;0.5 0.5 0.5;0.7 0.7 0.7]);
for i = 1:6
    load([datdir files{i}])
    fn = fieldnames(Avgs);
    if length(fn) == 3
        g = 1:3;
    else
        g = 1;
    end
    for j = g
        u = nanmean(Avgs.(fn{j}).x,2);
        z = hab(i,j)-0.04-linspace(0,0.03,35);
        zhc = z./hc(i);
        markx = u(1:5:end);
        marky = zhc(1:5:end);
        sp(i) = subplot(3,2,s(i));
        if i < 4
            plot(u,zhc,line{j},...
            'color',cl(j,:),...
            'linewidth',1.5), hold on
            plot(markx,marky,symb{j},...
                'linewidth',1.5,...
                'markerfacecolor',cl(j,:),...
                'markeredgecolor','k',...
                'markersize',6), hold on
        else
            plot(u,zhc,'-',...
            'color','k',...
            'linewidth',1.5), hold on
            plot(markx,marky,'^',...
                'linewidth',1.5,...
                'markerfacecolor','w',...
                'markeredgecolor','k',...
                'markersize',6)
        end
        grid on
    end
end
%global plot adjustments
set([sp(1) sp(2) sp(3)],'xlim',[-0.02 0.04],'gridlinestyle',':')
set([sp(4) sp(5) sp(6)],'xlim',[-0.015 0.005],'gridlinestyle',':')
set([sp(2) sp(3) sp(5) sp(6)],'xticklabel',[])
%HTA
set(sp(1),'ylim',[0 0.04],...
    'ytick',0:0.01:0.04,...
    'position',[0.13 0.1 0.35 0.25])
set(sp(2),'ylim',[0.29 0.35],...
    'ytick',0.29:0.02:0.35,...
    'position',[0.13 0.4 0.35 0.25])
set(sp(3),'ylim',[0.82 0.88],...
    'ytick',0.82:0.02:0.88,...
    'position',[0.13 0.7 0.35 0.25])
title(sp(3),'HTA')
ylabel(sp(2),'z/h_c')
xl = xlabel(sp(1),'$\overline{u} \quad (ms^{-1})$');
    set(xl,'Interpreter','latex','fontname','arial')
%VTA
set(sp(4),'ylim',[0.01 0.05],...
    'ytick',0.01:0.01:0.05,...
    'position',[0.6 0.1 0.35 0.25])
set(sp(5),'ylim',[0.58 0.625],...
    'ytick',0.58:0.015:0.63,...
    'position',[0.6 0.4 0.35 0.25])
set(sp(6),'ylim',[1.22 1.28],...
    'ytick',1.22:0.02:1.28,...
    'position',[0.6 0.7 0.35 0.25])
title(sp(6),'VTA')
xl = xlabel(sp(4),'$\overline{u} \quad (ms^{-1})$');
    set(xl,'Interpreter','latex','fontname','arial')
prettyfigures('text',13,'labels',14,'box',1,'tlength',[0.025 0.025])
export_fig([savefigdir 'VelProf_HTA_VTA'],'-pdf')

