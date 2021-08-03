%plot VPRO transects of NORMALIZED epsilon for all three deployments as subplots. 

clear
savefigs = 1;
dname = 'DPS2_depthTKE';
path = 'c:\Users\bkn5\Projects\Mekong_W2015\DataAnalysis\WavePowerFlux\DPS2Turbulence\';
D1.V1 = load([path 'Stat_VP1_140315.mat']);
D1.V2 = load([path 'Stat_VP2_140315.mat']);
D1.V3 = load([path 'Stat_VP3_140315.mat']);
mudflat = load([path 'WpfV5109_DPS2.mat']);
forest = load([path 'WpfAD5116_DPS2.mat']);
load('c:\Users\bkn5\Projects\Mekong_W2015\DataAnalysis\Spectra\DPS2\WaveStatsV5109_DPS2.mat')
azm = 360-wvstats.azm;
delx = 31; %m
theta = decimate((284-azm),11); %angle of waves relative to transect angle
h = decimate(forest.WPF.h,11)';
rho = decimate(forest.WPF.rho,11)';
wpf = abs((mudflat.WPF.Ffit-forest.WPF.Ffit)./(delx.*cosd(theta)))';

for i = 1:30
    D1.V1.E.beam1(i,:) = (D1.V1.STAT.TKE.Beam1(i,:).*h.*rho)./wpf;
    D1.V2.E.beam1(i,:) = (D1.V2.STAT.TKE.Beam1(i,:).*h.*rho)./wpf;
    D1.V3.E.beam1(i,:) = (D1.V3.STAT.TKE.Beam1(i,2:end).*h(2:end).*rho(2:end))./wpf(2:end);
    D1.V1.E.beam2(i,:) = (D1.V1.STAT.TKE.Beam2(i,:).*h.*rho)./wpf;
    D1.V2.E.beam2(i,:) = (D1.V2.STAT.TKE.Beam2(i,:).*h.*rho)./wpf;
    D1.V3.E.beam2(i,:) = (D1.V3.STAT.TKE.Beam2(i,2:end).*h(2:end).*rho(2:end))./wpf(2:end);
    D1.V1.E.beam3(i,:) = (D1.V1.STAT.TKE.Beam3(i,:).*h.*rho)./wpf;
    D1.V2.E.beam3(i,:) = (D1.V2.STAT.TKE.Beam3(i,:).*h.*rho)./wpf;
    D1.V3.E.beam3(i,:) = (D1.V3.STAT.TKE.Beam3(i,2:end).*h(2:end).*rho(2:end))./wpf(2:end);
    D1.V1.E.beam4(i,:) = (D1.V1.STAT.TKE.Beam4(i,:).*h.*rho)./wpf;
    D1.V2.E.beam4(i,:) = (D1.V2.STAT.TKE.Beam4(i,:).*h.*rho)./wpf;
    D1.V3.E.beam4(i,:) = (D1.V3.STAT.TKE.Beam4(i,2:end).*h(2:end).*rho(2:end))./wpf(2:end);
end
[n,m] = size(D1.V1.STAT.TKE.Beam1);
%calculate min and max for range plot
fid = fieldnames(D1);
for i = 1:length(fid)
    for ii = 1:n
        D1.(fid{i}).E.min(ii,:) = min(min(D1.(fid{i}).E.beam1(ii,:)),min(D1.(fid{i}).E.beam3(ii,:)));
        D1.(fid{i}).E.max(ii,:) = max(max(D1.(fid{i}).E.beam1(ii,:)),max(D1.(fid{i}).E.beam3(ii,:)));
    end  
end
f1 = figure;
set(f1,'PaperOrientation','portrait',...
    'position',[200 60   600 1000]);
set(gcf,'color','w','PaperPositionMode','auto')
ax = tight_subplot(3,1,0.025,[0.1 0.1],[0.125 0.5]);
colormap(paruly)
c = paruly(m);

axes(ax(3))
for i = 1:m
plot(log10((D1.V1.E.beam1(:,i)+D1.V1.E.beam3(:,i))./2),0.07-D1.V1.STAT.rangebins(1:30),...
    'Color',c(i,:),'Linewidth',1.5,'Marker','o',...
    'MarkerFaceColor',[1 1 1],'MarkerSize',5)
hold on
end
% plot(mean(D1.V1.E.beam3,2),D1.V1.STAT.rangebins(1:30),...
%     'Color','k','Linewidth',1.5,'Marker','^',...
%     'MarkerFaceColor',[1 1 1],'MarkerSize',10)
% plot(D1.V1.E.min,D1.V1.STAT.rangebins(1:30),...
%     'Color',[0.4 0.4 0.4],'LineStyle','--',...
%     'LineWidth',1.5)
% plot(D1.V1.E.max,D1.V1.STAT.rangebins(1:30),...
%     'Color',[0.4 0.4 0.4],'LineStyle','--',...
%     'LineWidth',1.5)
set(gca,'YDir','reverse')
% title('VP1 - 7cm','FontSize',20,'FontName','Arial')
xlabel(['\bf\it' '$Log_{10}\left( \hat{\epsilon}_{TKE}\right) $'],'interpreter','latex','FontSize',18)
grid on

axes(ax(2))
for i = 1:m
plot(log10((D1.V2.E.beam1(:,i)+D1.V2.E.beam3(:,i))./2),0.4-D1.V2.STAT.rangebins(1:30),...
    'Color',c(i,:),'Linewidth',1.5,'Marker','o',...
    'MarkerFaceColor',[1 1 1],'MarkerSize',5)
set(gca,'xticklabel',[])
hold on
end
% plot(mean(D1.V2.E.beam3,2),D1.V2.STAT.rangebins(1:30),...
%     'Color','k','Linewidth',1.5,'Marker','^',...
%     'MarkerFaceColor',[1 1 1],'MarkerSize',10)
% plot(D1.V2.E.min,D1.V2.STAT.rangebins(1:30),...
%     'Color',[0.4 0.4 0.4],'LineStyle','--',...
%     'LineWidth',1.5)
% plot(D1.V2.E.max,D1.V2.STAT.rangebins(1:30),...
%     'Color',[0.4 0.4 0.4],'LineStyle','--',...
%     'LineWidth',1.5)
set(gca,'YDir','reverse')
% xlabel('TKE Dissipation Rate (m^2/s^3)','FontSize',20,'FontName','Arial')
% title('VP2 - 40cm','FontSize',20,'FontName','Arial')
ylabel('Height Above Bed (m)','FontSize',20,'FontName','Arial')
cb = colorbar('eastoutside');ylabel(cb,'Minutes After Low Tide','FontSize',20)
caxis([0 m*10])
grid on

axes(ax(1))
for i = 1:m-1
plot(log10((D1.V3.E.beam1(:,i)+D1.V3.E.beam3(:,i))./2),0.8-D1.V3.STAT.rangebins(1:30),...
    'Color',c(i,:),'Linewidth',1.5,'Marker','o',...
    'MarkerFaceColor',[1 1 1],'MarkerSize',5)
set(gca,'xticklabel',[])
hold on
end

% plot(mean(D1.V3.E.beam3,2),D1.V3.STAT.rangebins(1:30),...
%     'Color','k','Linewidth',1.5,'Marker','^',...
%     'MarkerFaceColor',[1 1 1],'MarkerSize',10)
% plot(D1.V3.E.min,D1.V3.STAT.rangebins(1:30),...
%     'Color',[0.4 0.4 0.4],'LineStyle','--',...
%     'LineWidth',1.5)
% plot(D1.V3.E.max,D1.V3.STAT.rangebins(1:30),...
%     'Color',[0.4 0.4 0.4],'LineStyle','--',...
%     'LineWidth',1.5)
set(gca,'YDir','reverse')
% leg = legend({'Beam 1','Beam 3','Minimum','Maximum'},...
%     'Position',[0.78 0.13 0.12 0.1]);
% title('VP3 - 80cm','FontSize',20,'FontName','Arial')
grid on

%global plot adjustments
% set(ax(3),'YLim',[0.04 0.07],'YTick',0.04:0.005:0.07,...
%     'XLim',[-1 1],'XTick',-1:0.5:1,'GridLineStyle',':',...
%     'FontSize',20,'FontName','Arial','box','on')
set(ax,'XLim',[-2.5 1],'XTick',-2.5:0.5:1,'GridLineStyle',':',...
    'FontSize',20,'FontName','Arial','box','on')
set(ax,'LineWidth',1.5)
set(ax(1),'position',[0.17 0.7 0.48 0.26],'YLim',[0.73 0.76],'ytick',0.73:0.01:0.76,'ydir','normal')
set(ax(2),'position',[0.17 0.4 0.48 0.26],'YLim',[0.33 0.36],'ytick',0.33:0.01:0.36,'ydir','normal')
set(ax(3),'position',[0.17 0.1 0.48 0.26],'YLim',[0 0.03],'ytick',0:0.01:0.03,'ydir','normal')
set(cb,'position',[0.8 0.1 0.05 0.86],'linewidth',1.5)
% caxis([0 m*10])
if savefigs
    fname = [dname '_compare'];
    export_fig([path fname],'-pdf','-nocrop','-m1','-r900')
    disp(['Figure ' fname '.pdf saved'])
end

