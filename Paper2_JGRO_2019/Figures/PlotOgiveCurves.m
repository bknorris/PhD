%Plot cospectra, ogive curves from the CF method for estimating Reynolds
%Stress
clear,close all
load('d:\Projects\Mekong_W2015\DataAnalysis\Paper2\RS_estimates\RS_HTA1_vp3_sample.mat')
f1 = figure(1);
set(f1,'PaperOrientation','portrait',...
    'position',[400 100   800   600]);
set(gcf,'color','w','paperpositionmode','auto') 
%plot autospectra
sp(1) = subplot(311);
b(1) = loglog(ku,Su1,'color','k');hold on
b(2) = loglog(ku,Su2,'color',[0.6 0.6 0.6]);
l(1) = plot(repmat(kcu,10,1),linspace(10^-10,10^1,10),'color','k');
text(kcu+200,10^-2,'\bf\itf_{wc}')
lb = legend(b,{'Su_1u_1';'Su_2u_2'},'location','southwest');
%plot cospectrum
sp(2) = subplot(312);
p(1) = semilogx(ku,COuw,'color','k');hold on
l(2) = plot(repmat(kcu,10,1),linspace(-1,1,10),'color','k');
text(kcu+200,2E-3,'\bf\itf_{wc}')
%ogive curves
sp(3) = subplot(313);
dk = ku(4)-ku(3);
q(1) = semilogx(ku,cumtrapz(COuw),'color','k');hold on
q(2) = semilogx(ku,cumtrapz(COuwstar)/4.2,'color',[0.6 0.6 0.6]);
l(3) = plot(repmat(kcu,10,1),linspace(0,1,10),'color','k');
text(kcu+200,0.5,'\bf\itf_{wc}')
lq = legend(q,{'Integrated cospectrum';'Integrated model cospectrum'},'location','northwest');
%Global adjustments
set([b p q l],'linewidth',1.5)
set(sp(1),'position',[0.14 0.7 0.8 0.25],...
    'xticklabel',[],...
    'ylim',[10^-6 10^-1])
set(sp(2),'position',[0.14 0.41 0.8 0.25],...
    'xticklabel',[],...
    'ylim',[-1E-3 6.5E-3])
set(sp(3),'position',[0.14 0.12 0.8 0.25],...
    'ylim',[0 1])
set(sp,'xlim',[10 max(ku)])
%labels
txt = '$\int{Co_{u''w''}dk} (m/s)$';
xlabel(sp(3),'k (rad/m)')
ylabel(sp(3),txt,'interpreter','latex')
txt = 'Co_{u''w''}\n (m^2/s^2)/(rad/s) * 10^{-3}';
ylabel(sp(2),sprintf(txt))
txt = 'Autospectra\n (m^2/s^2)/(rad/s)';
ylabel(sp(1),sprintf(txt))
prettyfigures('text',12,'labels',14,'box',1)
set([lb lq],'box','off')
figdir = 'f:\GradSchool\DataAnalysis\Paper2\WorkingFigures\Spectra\';
export_fig([figdir 'COuw_ogive_HTA1'],'-pdf','-nocrop')

