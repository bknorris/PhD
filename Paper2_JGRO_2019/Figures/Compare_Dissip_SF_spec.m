%Plot dissipation rates generated from spectral means and from the
%structure function together to check for consistency. This script is for
%the thesis review.
clear, close all
maindir = 'e:\Mekong_W2015\DataAnalysis\Paper2\Turbulence\';
figdir = 'd:\GradSchool\DataAnalysis\Paper2\WorkingFigures\Turbulence\';
load('e:\Mekong_W2015\DataAnalysis\Paper2\Turbulence\Dissip_spectral_lcf_ucf.mat')
SFfiles = '8March2015_VelsTKE.mat';

f1 = figure(1);
set(f1,'PaperOrientation','portrait',...
    'position',[400 100   600   500],...
    'renderer','painters');
symb = {'^';'o';'s'};
c = csvread('d:\GradSchool\Design\Qualitative8_4.csv')./255;
days = 1;
sp = zeros(3,1);
load([maindir SFfiles])
plot(linspace(10^-5,10^-3,10),linspace(10^-5,10^-3,10),'-k','linewidth',1.5),hold on
hold on
for i = 1:3
    fn = fieldnames(Stat);
    dfn = fieldnames(TKE);
    bins = 13:17;
    Esp = TKE.(dfn{days}).(fn{i}).eww;Esp = nanmean(Esp(bins,:));
    Esf = (Stat.(fn{i}).z1.E+Stat.(fn{i}).z2.E)./2;Esf = nanmean(Esf(bins,:));
%     isimag = find(imag(Esp));
%     Esp(isimag) = NaN;
    fprintf('Esf is %0.2f%% of Esp\n',nanmean(Esf/Esp))
    sp(i) = plot(Esp,Esf,symb{i},...
        'color',c(i,:),...
        'linewidth',1.5,...
        'markersize',4,...
        'markerfacecolor',c(i,:));
end
leg = legend(sp,{'x = -10 cm';'x = 10 cm';'x = 20 cm'});
set(leg,'position',[0.23 0.82 0.05 0.05],...
    'box','off')
set(gca,'xscale','log',...
    'yscale','log')
axis equal
xlabel('\epsilon_{Sp} (W/kg)')
ylabel('\epsilon_{SF} (W/kg)')
title('HTA2 - March 8th, 2015')
prettyfigures('text',13,'labels',14,'box',1)
export_fig([figdir 'HTA2_dissip_compare'],'-pdf','-nocrop')