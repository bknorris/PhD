%Plot time-averaged velocity profiles for the HTA and VTA experiments. This
%is version 3 of this script
clear
datdir = 'f:\Mekong_W2015\DataAnalysis\Paper2\AveragedVelocities\SingleRotation\';
files = {'7March2015_Vave_sngrot.mat';'8March2015_Vave_sngrot.mat';'10March2015_Vave_sngrot.mat';...
    '14March2015a_bdadj_srot.mat'};
savefigdir = 'g:\GradSchool\DataAnalysis\Paper2\WorkingFigures\VerticalProfiles\';

%plot routine
f1 = figure;
set(f1,'PaperOrientation','portrait',...
    'position',[400 100   800  800]);
set(gcf,'color','w','paperpositionmode','auto')
s = [7 7 7;4 4 4;1 1 1;9 6 3]; %subplot order
sp = zeros(9,1);
hab = [0.065 0.065 0.061;
    0.242 0.240 0.240;
    0.550 0.550 0.550;
    0.07  0.416 0.806];
hc = [0.58 0.58 0.58 0.6 0.6 0.6];
symb = {'o';'d';'p'};
line = {'-';'--';'-.'};
cl = flipud([0.1 0.1 0.1;0.5 0.5 0.5;0.7 0.7 0.7]);
for i = 1:4
    load([datdir files{i}])
    fn = fieldnames(Avgs);
    for j = 1:3
        u = nanmean(Avgs.(fn{j}).x,2);
        w = nanmean(((Avgs.(fn{j}).z1+Avgs.(fn{j}).z2)./2),2);
        if i > 3
            sb = '^';cc = [0 0 0];
        else
            sb = symb{j};cc = cl(j,:);
        end
        if i == 1 || i == 4
            u(25:35) = NaN;w(25:35) = NaN;
            u(25:35) = NaN;w(25:35) = NaN;
        end
        z = hab(i,j)-0.04-linspace(0,0.03,35);
        zhc = z./hc(i);
        markxu = u(1:5:end);
        markxw = w(1:8:end);
        markyu = zhc(1:5:end);
        markyw = zhc(1:8:end);
        sp(s(i,j)) = subplot(3,3,s(i,j));
        plot(u,zhc,'-',...
            'color',cc,...
            'linewidth',1.5), hold on
        plot(markxu,markyu,sb,...
            'linewidth',1.5,...
            'markerfacecolor',cc,...
            'markeredgecolor','k',...
            'markersize',5)
        plot(w,zhc,'-.',...
            'color',cc,...
            'linewidth',1.5)
        plot(markxw,markyw,sb,...
            'linewidth',1.5,...
            'markerfacecolor',cc,...
            'markeredgecolor','k',...
            'markersize',5)
        grid on
    end
end

%plot vorticity in only HTA experiments
s = [8 5 2]; %subplot order
for i = 1:3
    load([datdir files{i}])
    %calc omega_y between VP2 and VP3
    w1 = (Avgs.vpro2.z1+Avgs.vpro2.z2)./2;
    w2 = (Avgs.vpro3.z1+Avgs.vpro3.z2)./2;
    x1 = Avgs.vpro2.x;x2 = Avgs.vpro3.x;
    if length(w1) < length(w2)
        w2 = w2(:,1:length(w1));
        x2 = x2(:,1:length(x1));
    end
    dw = nanmean(w1-w2,2);
    if i == 1
        du = nanmean((x1+x2)./2,2);
    else
        du = nanmean(x1,2);
    end
    omgy = du-dw;
    z = hab(i,j)-0.04-linspace(0,0.03,35);
    zhc = z./hc(i);
    markxu = omgy(1:5:end);
    markyu = zhc(1:5:end);
    sp(s(i)) = subplot(3,3,s(i));
    plot(omgy,zhc,'-',...
        'color','k',...
        'linewidth',1.5), hold on
    plot(markxu,markyu,'o',...
        'linewidth',1.5,...
        'markerfacecolor','k',...
        'markeredgecolor','k',...
        'markersize',4)
    grid on
end
%global plot adjustments
set(sp,'gridlinestyle',':')
set(sp(1),'xlim',[-0.02 0.04],'xtick',-0.02:0.03:0.04)
set(sp(2),'xlim',[0.02 0.04],'xtick',0.02:0.01:0.04)
set(sp(3),'xlim',[-0.02 0],'xtick',-0.02:0.01:0)
set(sp(4),'xlim',[-0.02 0.02],'xtick',-0.02:0.02:0.02)
set(sp(5),'xlim',[0.005 0.015],'xtick',0.005:0.005:0.015)
set(sp(6),'xlim',[-0.01 0.005],'xtick',-0.01:0.0075:0.005)
set(sp(7),'xlim',[-0.02 0.02],'xtick',-0.02:0.02:0.02)
set(sp(8),'xlim',[-0.01 0.01],'xtick',-0.01:0.01:0.01)
set(sp(9),'xlim',[-0.005 0.015],'xtick',-0.005:0.01:0.015)
%HTA
%row 3
set(sp(7),'ylim',[0 0.04],...
    'ytick',0:0.01:0.04,...
    'position',[0.12 0.1 0.21 0.22])
set(sp(8),'ylim',[0 0.04],...
    'ytick',0:0.01:0.04,...
    'position',[0.44 0.1 0.21 0.22])
set(sp(9),'ylim',[0.01 0.05],...
    'ytick',0.01:0.01:0.05,...
    'position',[0.76 0.1 0.21 0.22])
%row 2
set(sp(4),'ylim',[0.29 0.35],...
    'ytick',0.29:0.02:0.35,...
    'position',[0.12 0.42 0.21 0.22])
set(sp(5),'ylim',[0.29 0.35],...
    'ytick',0.29:0.02:0.35,...
    'position',[0.44 0.42 0.21 0.22])
set(sp(6),'ylim',[0.58 0.625],...
    'ytick',0.58:0.015:0.63,...
    'position',[0.76 0.42 0.21 0.22])
%row 1
set(sp(1),'ylim',[0.82 0.88],...
    'ytick',0.82:0.02:0.88,...
    'position',[0.12 0.74 0.21 0.22])
set(sp(2),'ylim',[0.82 0.88],...
    'ytick',0.82:0.02:0.88,...
    'position',[0.44 0.74 0.21 0.22])
set(sp(3),'ylim',[1.22 1.28],...
    'ytick',1.22:0.02:1.28,...
    'position',[0.76 0.74 0.21 0.22])
title(sp(1),'(a)'),ylabel(sp(1),'\itz/h_c')
title(sp(2),'(d)')
title(sp(3),'(g)')
title(sp(4),'(b)'),ylabel(sp(4),'\itz/h_c')
title(sp(5),'(e)')
title(sp(6),'(h)')
title(sp(7),'(c)'),ylabel(sp(7),'\itz/h_c'),xlabel(sp(7),'u or w (ms^{-1})')
title(sp(8),'(f)'),xlabel(sp(8),'\omega_y (s^{-1})')
title(sp(9),'(i)'),xlabel(sp(9),'u or w (ms^{-1})')
prettyfigures('text',12,'labels',13,'box',1,'tlength',[0.025 0.025])
% export_fig([savefigdir 'VelProf_HTA_VTA_v5'],'-pdf')

