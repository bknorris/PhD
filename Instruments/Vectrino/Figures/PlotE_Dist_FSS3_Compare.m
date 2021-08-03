%plot VPRO transects of epsilon for all three deployments as subplots. ONR
%Talk.

clear

FSS1.V1 = load('Stat_VP1_070315.mat');
FSS1.V2 = load('Stat_VP2_070315.mat');
FSS1.V3 = load('Stat_VP3_070315.mat');

FSS2.V1 = load('Stat_VP1_080315.mat');
FSS2.V2 = load('Stat_VP2_080315.mat');
FSS2.V3 = load('Stat_VP3_080315.mat');

FSS3.V1 = load('Stat_VP1_100315.mat');
FSS3.V2 = load('Stat_VP2_100315.mat');
FSS3.V3 = load('Stat_VP3_100315.mat');

dname = 'FSS3';
savefigs = 1;
savefigdir = 'c:\Users\bkn5\Desktop\';
%distances are based on photographs of the FSS3 setup:
%VP1 is 10cm from the pneumatophores
%VP2 is 10cm from the pneumatophores and 20cm from VP1
%VP3 is 20cm from the pneumatophores and 30cm from VP1
%VP3 is closest to N, VP1 is closest to S
%Vectrino beams are 2.5cm apart from the center transducer
%define variables:
m = length(FSS1.V1.STAT.TKE.Time)-1;
intv = 10; %averaging interval
bin15 = zeros(m,7);
dist = [-12.5 -7.5 0 7.5 12.5 17.5 22.5];
for i = 1:m; %time
    FSS1.bin15(i,:) = [log10(FSS1.V1.STAT.TKE.Beam1(5,i))  log10(FSS1.V1.STAT.TKE.Beam4(5,i))...
        NaN log10(FSS1.V2.STAT.TKE.Beam1(5,i)) log10(FSS1.V2.STAT.TKE.Beam4(5,i))...
        log10(FSS1.V3.STAT.TKE.Beam1(5,i)) log10(FSS1.V3.STAT.TKE.Beam4(5,i))];
    FSS1.bin15(i,3) = ((FSS1.bin15(i,2) + FSS1.bin15(i,4))/2);
end

m = length(FSS2.V1.STAT.TKE.Time)-1;
for i = 1:m; %time
    FSS2.bin15(i,:) = [log10(FSS2.V1.STAT.TKE.Beam1(15,i))  log10(FSS2.V1.STAT.TKE.Beam4(15,i))...
        NaN log10(FSS2.V2.STAT.TKE.Beam1(15,i)) log10(FSS2.V2.STAT.TKE.Beam4(15,i))...
        log10(FSS2.V3.STAT.TKE.Beam1(15,i)) log10(FSS2.V3.STAT.TKE.Beam4(15,i))];
    FSS2.bin15(i,3) = ((FSS2.bin15(i,2) + FSS2.bin15(i,4))/2);
end

m = length(FSS3.V1.STAT.TKE.Time)-1;
for i = 1:m; %time
    FSS3.bin15(i,:) = [log10(FSS3.V1.STAT.TKE.Beam1(15,i))  log10(FSS3.V1.STAT.TKE.Beam4(15,i))...
        NaN log10(FSS3.V2.STAT.TKE.Beam1(15,i)) log10(FSS3.V2.STAT.TKE.Beam4(15,i))...
        log10(FSS3.V3.STAT.TKE.Beam1(15,i)) log10(FSS3.V3.STAT.TKE.Beam4(15,i))];
    FSS3.bin15(i,3) = ((FSS3.bin15(i,2) + FSS3.bin15(i,4))/2);
end


t = max(FSS2.V1.STAT.TKE.Time);

f1 = figure;
set(f1,'PaperOrientation','portrait',...
    'position',[500 60   800 1050]);
set(gcf,'color','w','PaperPositionMode','auto')
ax = tight_subplot(3,1,0.025,[0.1 0.1],[0.125 0.125]);
colormap(linspecer)

axes(ax(1))
m = length(FSS3.V1.STAT.TKE.Time)-1;
c = linspecer(18,'sequential');
for i = 1:m
plot(dist,FSS3.bin15(i,:),'-+','Color',c(i,:),'linewidth',1.5)
hold on
end
grid on
set(ax(1),'XTickLabel',[])
text(0.055,0.95,'\bf\itVP1','units','normalized','FontSize',14)
text(0.6,0.95,'\bf\itVP2','units','normalized','FontSize',14)
text(0.9,0.95,'\bf\itVP3','units','normalized','FontSize',14)
line([0 0],[-10 0],'Color',[0 0 0],'linewidth',2)
% title('\bf\itFSS3.1 Vectrino Cross-Section 60mm above bed','FontSize',14)
cb = colorbar('northoutside');
caxis([0 t])
set(get(cb,'title'),'string','\bf\itMinutes Elapsed','FontSize',14);

axes(ax(2))
m = length(FSS2.V1.STAT.TKE.Time)-1;
c = linspecer(m,'sequential');
for i = 1:m
plot(dist,FSS2.bin15(i,:),'-+','Color',c(i,:),'linewidth',1.5)
hold on
end
grid on
set(ax(2),'XTickLabel',[])
text(0.055,0.95,'\bf\itVP1','units','normalized','FontSize',14)
text(0.6,0.95,'\bf\itVP2','units','normalized','FontSize',14)
text(0.9,0.95,'\bf\itVP3','units','normalized','FontSize',14)
ylabel('\bf\itlog_1_0(\epsilon)','FontSize',14)
line([0 0],[-10 0],'Color',[0 0 0],'linewidth',2)

axes(ax(3))
m = length(FSS1.V1.STAT.TKE.Time)-1;
c = linspecer(m,'sequential');
for i = 1:m
plot(dist,FSS1.bin15(i,:),'-+','Color',c(i,:),'linewidth',1.5)
hold on
end
grid on
text(0.055,0.95,'\bf\itVP1','units','normalized','FontSize',14)
text(0.6,0.95,'\bf\itVP2','units','normalized','FontSize',14)
text(0.9,0.95,'\bf\itVP3','units','normalized','FontSize',14)
xlabel('\bf\itDistance Along Cross-Shore Transect (cm)','FontSize',14)
line([0 0],[-10 0],'Color',[0 0 0],'linewidth',2)

%global plot adjustments
set(ax,'YLim',[-6 -1],'YTick',-6:1:-1,...
    'XLim',[-13 23],'XTick',-13:3:23,'GridLineStyle',':',...
    'FontSize',14,'FontName','Arial')
set(ax(1),'position',[0.1 0.63 0.82 0.25])
set(ax(2),'position',[0.1 0.355 0.82 0.25])
set(ax(3),'position',[0.1 0.08 0.82 0.25])
if savefigs
    fpath = savefigdir;fname = [dname '_compare'];
    export_fig([fpath fname],'-png','-m1','-r900','-opengl')
    disp(['Figure ' fname '.pdf saved'])
end

