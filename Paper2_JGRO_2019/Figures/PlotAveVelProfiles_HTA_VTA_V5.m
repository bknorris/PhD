%Plot time-averaged velocity profiles for the HTA and VTA experiments. This
%is version 5 of this script
clear
datdir = 'D:\Projects\Mekong_W2015\DataAnalysis\Paper2\AveragedVelocities\Two\';
files = {'7March2015_Vave.mat';'8March2015_Vave.mat';'10March2015_Vave.mat';...
    '14March2015a_Vave.mat'};
savefigdir = 'g:\GradSchool\DataAnalysis\Paper2\WorkingFigures\VerticalProfiles\';
cmap = csvread('g:\GradSchool\Design\Qualitative8_4.csv');
cmap = flipud(cmap)./255;

%plot routine
f1 = figure;
set(f1,'PaperOrientation','portrait',...
    'position',[400 100   600  700]);
set(gcf,'color','w','paperpositionmode','auto')
s = [5 5 5;3 3 3;1 1 1;6 4 2]; %subplot order
sp = zeros(6,1);
hab = [0.066 0.066 0.066;
    0.248 0.248 0.248;
    0.550 0.550 0.550;
    0.07  0.416 0.806];
hc = [0.58 0.58 0.58 0.6 0.6 0.6];
symb = {'o';'d';'p'};
line = {'-';'--';'-.'};
for i = 1:4
    load([datdir files{i}])
    fn = fieldnames(Avgs);
    for j = 1:3
        if i == 4
            bid = 14;
        else
            bid = 1;
        end
        u = nanmean(Avgs.(fn{j}).x(:,bid:end),2);
        w = nanmean(((Avgs.(fn{j}).z1(:,bid:end)+Avgs.(fn{j}).z2(:,bid:end))./2),2);
        if i < 4
            cc = cmap(j,:);
        else
            cc = cmap(4,:);
        end
        if i == 1 || i == 4
            u(25:35) = NaN;w(25:35) = NaN;
            u(25:35) = NaN;w(25:35) = NaN;
        end
        z = hab(i,j)-0.04-linspace(0,0.03,35);
        zhc = z./hc(i);
%         markxu = u(1:5:end);
%         markxw = w(1:8:end);
%         markyu = zhc(1:5:end);
%         markyw = zhc(1:8:end);
        sp(s(i,j)) = subplot(3,2,s(i,j));
        plot(u,zhc,'-',...
            'color',cc,...
            'linewidth',3), hold on
%         plot(markxu,markyu,sb,...
%             'linewidth',1.5,...
%             'markerfacecolor',cc,...
%             'markeredgecolor','k',...
%             'markersize',5)
        plot(w,zhc,'-.',...
            'color',cc,...
            'linewidth',3)
%         plot(markxw,markyw,sb,...
%             'linewidth',1.5,...
%             'markerfacecolor',cc,...
%             'markeredgecolor','k',...
%             'markersize',5)
        %plot horizontal line in subplot 323
        if i == 2
            plot(linspace(-0.02,0.04,10),ones(10,1)*0.33,'--k','linewidth',1.5)
        end
        grid on
    end
end
%global plot adjustments
set(sp,'gridlinestyle',':')
set([sp(1) sp(3) sp(5)],'xlim',[-0.02 0.04],'xtick',-0.02:0.02:0.04)
set([sp(2) sp(4) sp(6)],'xlim',[-0.02 0.015],'xtick',-0.02:0.01:0.01)
set(sp(3),'xlim',[-0.02 0.04],'xtick',-0.02:0.02:0.04)
%HTA
set(sp(1),'ylim',[0.82 0.88],'ytick',0.82:0.02:0.88,'position',[0.15 0.7 0.32 0.21])
set(sp(3),'ylim',[0.3 0.36],'ytick',0.3:0.02:0.36,'position',[0.15 0.4 0.32 0.21])
set(sp(5),'ylim',[0 0.05],'ytick',0:0.02:0.6,'position',[0.15 0.1 0.32 0.21])
%VTA
set(sp(2),'ylim',[1.24 1.28],'ytick',1.24:0.02:1.28,'position',[0.58 0.7 0.32 0.21])
set(sp(4),'ylim',[0.59 0.63],'ytick',0.58:0.02:0.64,'position',[0.58 0.4 0.32 0.21])
set(sp(6),'ylim',[0 0.05],'ytick',0:0.02:0.6,'position',[0.58 0.1 0.32 0.21])
title(sp(1),'(a)'),ylabel(sp(1),'\itz/h_c')
title(sp(3),'(b)'),ylabel(sp(3),'\itz/h_c')
title(sp(5),'(c)'),ylabel(sp(5),'\itz/h_c'),xlabel(sp(5),'{u or w} (m/s)')
title(sp(2),'(d)')
title(sp(4),'(e)')
title(sp(6),'(f)'),xlabel(sp(6),'{u or w} (m/s)')
% prettyfigures('text',12,'labels',13,'box',1,'tlength',[0.025 0.025])
% export_fig([savefigdir 'VelProf_HTA_VTA_v6'],'-pdf')

