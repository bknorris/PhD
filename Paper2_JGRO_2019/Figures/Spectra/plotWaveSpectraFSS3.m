%make spectra plot for the Student Conference talk
cd('/Volumes/MEKONG1/Mekong_W2015/DataAnalysis/WavePowerFlux/FSS3.1')
fss1p = load('V5108psd.mat');
fss1f = load('V5108freq.mat');
cd ../FSS3.2/
fss2p = load('HR3psd.mat');
fss2f = load('HR3freq.mat');
cd ../FSS3.3/
fss3p = load('HR3psd.mat');
fss3f = load('HR3freq.mat');
cd ../

fss1p.pdens(1:200) = 0;
fss2p.pdens(1:200) = 0;
fss3p.pdens(1:200) = 0;

figure
set(gcf,'color','w','PaperPositionMode','auto'), hold on
p(3) = line_fewer_markers(fss3f.freq,fss3p.pdens,50,'x-k','Spacing','Curve',...
    'MarkerSize', 7, 'MarkerFaceColor', 'k','MarkerEdgeColor','k',...
    'LineWidth',1.5,'LegendLine',1);
p(2) = line_fewer_markers(fss2f.freq,fss2p.pdens,50,'^-r','Spacing','Curve',...
    'MarkerSize', 6, 'MarkerFaceColor', 'r','MarkerEdgeColor','r',...
    'LineWidth',1.5,'LegendLine',1);
p(1) = line_fewer_markers(fss1f.freq,fss1p.pdens,50,'.-b','Spacing','Curve',...
    'MarkerSize', 15, 'MarkerFaceColor', 'b','MarkerEdgeColor','b',...
    'LineWidth',1.5,'LegendLine',1);
xlabel('\bf\itFrequency (Hz)','FontSize',14)
ylabel('\bf\itPower (m^2/Hz)','FontSize',14)
ax1 = gca;
set(ax1,'Xlim',[0.05 1],'FontSize',12,'LineWidth',1.5)
leg = legend(p,'\bf\it07-03-2015','\bf\it08-03-2015','\bf\it10-03-2015',...
    'location','northeast');
legend boxoff
box on
fpwd = '/Volumes/MEKONG1/Mekong_W2015/Figures/Spectra/';
export_fig([fpwd 'FSS3spectra'],'-pdf','-nocrop','-m1','-r900')

