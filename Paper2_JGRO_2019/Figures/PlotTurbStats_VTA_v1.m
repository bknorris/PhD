%Plot averaged Reynolds Stresses and eddy viscosities
%for the VTA experiment (2 x 3 panel subplots). This is version 1 of this
%script.
clear, close all
cmap = csvread('f:\GradSchool\Design\Qualitative8_4.csv');
cmap = flipud(cmap)./255;

%plot routine
f1 = figure;
set(f1,'PaperOrientation','portrait',...
    'position',[400 100   600  700]);
set(gcf,'color','w','paperpositionmode','auto')
sp = zeros(6,1);
hab = [0.07  0.416 0.806];

hc = 0.6;
savefigdir = 'f:\GradSchool\DataAnalysis\Paper2\WorkingFigures\VerticalProfiles\';
%% Load the Reynolds stresses
load('D:\Projects\Mekong_W2015\DataAnalysis\Paper2\RS_estimates\RStress_10min_v4test.mat')
s = [5 3 1];
fn = fieldnames(RS.day4);
for j = 1:3
    tmp = RS.day4.(fn{j}).uw;
    tmp = tmp(15:40,:);
    id = find(tmp> 0.25 | tmp < -0.25);
    tmp(id) = NaN;
    uw = nanmean(tmp);
    z = hab(j)-0.04-linspace(0,0.03,35);
    zhc = z./hc;
    sp(s(j)) = subplot(3,2,s(j));
    plot(uw,zhc,'-',...
        'color',cmap(4,:),...
        'linewidth',3), hold on
    grid on
end
%% Load the eddy viscosity
load('D:\Projects\Mekong_W2015\DataAnalysis\Paper2\EddyViscosity\U_z_gradient_v4.mat')
s = [6 4 2];
for i = 1:4
    fn = fieldnames(RS.day4);
    for j = 1:3
        uw = RS.day4.(fn{j}).uw;
        time1 = RS.day4.(fn{j}).time;
        time2 = slope.day4.(fn{j}).time';
        dudz = repmat(slope.day4.(fn{j}).dudz,1,35);
        %eddy viscosity is Re stress/velocity gradient
        nue = -1*nanmean(uw)./dudz;
        nue(isinf(nue)) = NaN;
        z = hab(j)-0.04-linspace(0,0.03,35);
        zhc = z./hc;
            sp(s(j)) = subplot(3,2,s(j));
            plot(nue,zhc,'-',...
                'color',cmap(4,:),...
                'linewidth',3), hold on
            grid on
    end
end
%% Global plot adjustments
set(sp,'gridlinestyle',':')
set([sp(1) sp(3) sp(5)],'xlim',[-2E-4 6E-4],'xtick',-2E-4:2E-4:6E-4)
set([sp(2) sp(4) sp(6)],'xlim',[-2E-4 4E-4],'xtick',-2E-4:2E-4:4E-4)
%HTA
set(sp(1),'ylim',[1.24 1.28],'ytick',1.24:0.02:1.28,'position',[0.15 0.7 0.32 0.21])
set(sp(3),'ylim',[0.59 0.63],'ytick',0.58:0.02:0.64,'position',[0.15 0.4 0.32 0.21])
set(sp(5),'ylim',[0 0.05],'ytick',0:0.02:0.6,'position',[0.15 0.1 0.32 0.21])
%VTA
set(sp(2),'ylim',[1.24 1.28],'ytick',1.24:0.02:1.28,'position',[0.58 0.7 0.32 0.21])
set(sp(4),'ylim',[0.59 0.63],'ytick',0.58:0.02:0.64,'position',[0.58 0.4 0.32 0.21])
set(sp(6),'ylim',[0 0.05],'ytick',0:0.02:0.6,'position',[0.58 0.1 0.32 0.21])
title(sp(1),'(a)'),ylabel(sp(1),'\itz/h_c')
title(sp(3),'(b)'),ylabel(sp(3),'\itz/h_c')
title(sp(5),'(c)'),ylabel(sp(5),'\itz/h_c'),xlabel(sp(5),'$-\overline{u''w''} \quad (m^2s^{-2})$','Interpreter','latex')
title(sp(2),'(d)')
title(sp(4),'(e)')
title(sp(6),'(f)'),xlabel(sp(6),'$\overline{\nu_t} \quad (m^2s^{-1})$','Interpreter','latex')
prettyfigures('text',12,'labels',13,'box',1,'tlength',[0.025 0.025])
export_fig([savefigdir 'TurbStats_VTA_v1'],'-pdf')