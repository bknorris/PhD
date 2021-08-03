%plot VPRO transects of NORMALIZED epsilon for all three deployments as subplots. 

clear
pwd = 'c:/users/bkn5/Projects/Mekong_W2015/DataAnalysis/WavePowerFlux/FSS3Turbulence/';
FSS1.V1 = load([pwd 'Stat_VP1_070315.mat']);
FSS1.V2 = load([pwd 'Stat_VP2_070315.mat']);
FSS1.V3 = load([pwd 'Stat_VP3_070315.mat']);
mudflat = load([pwd 'WpfV5108_FSS31.mat']);
forest = load([pwd 'WpfAD5116_FSS31.mat']);
delx = 65; %m
theta = (296-350); %angle of V5108 to AD5116 relative to transect angle
h = decimate(forest.WPF.h,10)';
rho = decimate(forest.WPF.rho,10)';
wpf = abs((mudflat.WPF.Ffit-forest.WPF.Ffit)./(delx.*cosd(theta)))';

for i = 1:30
    FSS1.V1.E.beam1(i,:) = (FSS1.V1.STAT.TKE.Beam1(i,8:17).*h.*rho)./wpf;
    FSS1.V2.E.beam1(i,:) = (FSS1.V2.STAT.TKE.Beam1(i,8:17).*h.*rho)./wpf;
    FSS1.V3.E.beam1(i,:) = (FSS1.V3.STAT.TKE.Beam1(i,8:17).*h.*rho)./wpf;
    FSS1.V1.E.beam3(i,:) = (FSS1.V1.STAT.TKE.Beam3(i,8:17).*h.*rho)./wpf;
    FSS1.V2.E.beam3(i,:) = (FSS1.V2.STAT.TKE.Beam3(i,8:17).*h.*rho)./wpf;
    FSS1.V3.E.beam3(i,:) = (FSS1.V3.STAT.TKE.Beam3(i,8:17).*h.*rho)./wpf;
end
pwd = 'c:/users/bkn5/Projects/Mekong_W2015/DataAnalysis/WavePowerFlux/FSS3Turbulence/';
FSS2.V1 = load([pwd 'Stat_VP1_080315.mat']);
FSS2.V2 = load([pwd 'Stat_VP2_080315.mat']);
FSS2.V3 = load([pwd 'Stat_VP3_080315.mat']);
mudflat = load([pwd 'WpfADHR3_FSS32.mat']);
forest = load([pwd 'WpfAD5116_FSS32.mat']);
delx = 56; %m
theta = 0; 
h = decimate(forest.WPF.h,10)';
rho = decimate(forest.WPF.rho,10)';
wpf = abs((mudflat.WPF.Ffit-forest.WPF.Ffit)./(delx.*cosd(theta)))';
for i = 1:30
    FSS2.V1.E.beam1(i,:) = (FSS2.V1.STAT.TKE.Beam1(i,1:23).*h.*rho)./wpf;
    FSS2.V2.E.beam1(i,:) = (FSS2.V2.STAT.TKE.Beam1(i,1:23).*h.*rho)./wpf;
    FSS2.V3.E.beam1(i,:) = (FSS2.V3.STAT.TKE.Beam1(i,1:23).*h.*rho)./wpf;
    FSS2.V1.E.beam3(i,:) = (FSS2.V1.STAT.TKE.Beam3(i,1:23).*h.*rho)./wpf;
    FSS2.V2.E.beam3(i,:) = (FSS2.V2.STAT.TKE.Beam3(i,1:23).*h.*rho)./wpf;
    FSS2.V3.E.beam3(i,:) = (FSS2.V3.STAT.TKE.Beam3(i,1:23).*h.*rho)./wpf;
end
pwd = 'c:/users/bkn5/Projects/Mekong_W2015/DataAnalysis/WavePowerFlux/FSS3Turbulence/';
FSS3.V1 = load([pwd 'Stat_VP1_100315.mat']);
FSS3.V2 = load([pwd 'Stat_VP2_100315.mat']);
FSS3.V3 = load([pwd 'Stat_VP3_100315.mat']);
mudflat = load([pwd 'WpfHR3_FSS33.mat']);
forest = load([pwd 'WpfAD5116_FSS33.mat']);
delx = 56; %m
theta = 0; 
h = decimate(forest.WPF.h,11)';
rho = decimate(forest.WPF.rho,11)';
wpf = abs((mudflat.WPF.Ffit-forest.WPF.Ffit)./(delx.*cosd(theta)))';
for i = 1:30
    FSS3.V1.E.beam1(i,:) = (FSS3.V1.STAT.TKE.Beam1(i,:).*h.*rho)./wpf;
    FSS3.V2.E.beam1(i,:) = (FSS3.V2.STAT.TKE.Beam1(i,:).*h.*rho)./wpf;
    FSS3.V3.E.beam1(i,:) = (FSS3.V3.STAT.TKE.Beam1(i,:).*h.*rho)./wpf;
    FSS3.V1.E.beam3(i,:) = (FSS3.V1.STAT.TKE.Beam3(i,:).*h.*rho)./wpf;
    FSS3.V2.E.beam3(i,:) = (FSS3.V2.STAT.TKE.Beam3(i,:).*h.*rho)./wpf;
    FSS3.V3.E.beam3(i,:) = (FSS3.V3.STAT.TKE.Beam3(i,:).*h.*rho)./wpf;
end

dname = 'FSS3_longestTS';
savefigs = 0;
savefigdir = 'c:/users/bkn5/Projects/Mekong_W2015/Figures/Turbulence/';
%create time vector to orient the color spacing of the lines according to
%the inundation time of the offshore vector (V5108)
t0 = datenum(2015,03,07,13,22,30); %time water first submerges the mudflat Vector on Day 1 (beginning of time vector)
t1 = datenum(FSS1.V1.STAT.DepStop(1:end-4),'dd-mm-yyyy HH:MM:SS'); %final time of measurement (end of time vector)
vpt0 = datenum(FSS1.V1.STAT.DepStart(1:end-4),'dd-mm-yyyy HH:MM:SS'); %get time of vecpro start from Stat file
step = datenum(0,0,0,0,10,0); %timestep of 10mins
time = t0:step:t1;
FSS1.time = 1:length(time);
[~,FSS1.id] = min(abs(time-vpt0));
t0 = datenum(2015,03,08,13,42,30); %time water first submerges the mudflat Vector on Day 2 (beginning of time vector)
t1 = datenum(FSS2.V1.STAT.DepStop(1:end-4),'dd-mm-yyyy HH:MM:SS'); %final time of measurement (end of time vector)
vpt0 = datenum(FSS2.V1.STAT.DepStart(1:end-4),'dd-mm-yyyy HH:MM:SS');
time = t0:step:t1;
FSS2.time = 1:length(time);
[~,FSS2.id] = min(abs(time-vpt0));
t0 = datenum(2015,03,10,14,07,00); %time water first submerges the mudflat Vector on Day 3 (beginning of time vector)
t1 = datenum(FSS3.V1.STAT.DepStop(1:end-4),'dd-mm-yyyy HH:MM:SS'); %final time of measurement (end of time vector)
vpt0 = datenum(2015,03,10,15,20,00);
time = t0:step:t1;
FSS3.time = 1:length(time);
[~,FSS3.id] = min(abs(time-vpt0));
%distances are based on photographs of the FSS3 setup:
%VP1 is 10cm from the pneumatophores
%VP2 is 10cm from the pneumatophores and 20cm from VP1
%VP3 is 20cm from the pneumatophores and 30cm from VP1
%VP3 is closest to N, VP1 is closest to S
%Vectrino beams are 2.5cm apart from the center transducer
%define variables:
m = length(FSS1.V1.STAT.TKE.Time)-11;
intv = 10; %averaging interval
bin15 = zeros(m,10);
dist = [-12.5 -7.5 0 7.5 12.5 17.5 22.5];
for i = 1:m; %time
    FSS1.bin15(i,:) = [FSS1.V1.E.beam1(7,i)  FSS1.V1.E.beam3(7,i)...
        NaN FSS1.V2.E.beam1(7,i) FSS1.V2.E.beam3(7,i)...
        FSS1.V3.E.beam1(7,i) FSS1.V3.E.beam3(7,i)];
    FSS1.bin15(i,3) = ((FSS1.bin15(i,2) + FSS1.bin15(i,4))/2);
end

m = length(FSS2.V1.STAT.TKE.Time)-6;
for i = 1:m; %time
    FSS2.bin15(i,:) = [FSS2.V1.E.beam1(15,i)  FSS2.V1.E.beam3(15,i)...
        NaN FSS2.V2.E.beam1(15,i)  FSS2.V2.E.beam3(15,i)...
        FSS2.V3.E.beam1(15,i)  FSS2.V3.E.beam3(15,i)];
    FSS2.bin15(i,3) = ((FSS2.bin15(i,2) + FSS2.bin15(i,4))/2);
end

m = length(FSS3.V1.STAT.TKE.Time)-1;
for i = 1:m; %time
    FSS3.bin15(i,:) = [FSS3.V1.E.beam1(15,i)  FSS3.V1.E.beam3(15,i)...
        NaN FSS3.V2.E.beam1(15,i)  FSS3.V2.E.beam3(15,i)...
        FSS3.V3.E.beam1(15,i)  FSS3.V3.E.beam3(15,i)];
    FSS3.bin15(i,3) = ((FSS3.bin15(i,2) + FSS3.bin15(i,4))/2);
end

f1 = figure;
set(f1,'PaperOrientation','portrait',...
    'position',[500 60   800 1050]);
set(gcf,'color','w','PaperPositionMode','auto')
ax = tight_subplot(3,1,0.025,[0.1 0.1],[0.125 0.125]);
colormap(paruly)

axes(ax(1))
m = length(FSS2.time)-11;
n = length(FSS3.bin15);
c = paruly(m);
t = FSS3.id:m;
for i = 1:n
plot(dist,log10(FSS3.bin15(i,:)),'Color',c(t(i),:),'linewidth',2,...
    'Marker','o','MarkerFaceColor',[1 1 1],'MarkerSize',6)
hold on
end
grid on
set(ax(1),'XTickLabel',[])
text(0.055,0.92,'\bf\itVP1','units','normalized','FontSize',20)
text(0.6,0.92,'\bf\itVP2','units','normalized','FontSize',20)
text(0.9,0.92,'\bf\itVP3','units','normalized','FontSize',20)
line([0 0],[-2 4],'Color',[0 0 0],'linewidth',2)
% title('\bf\itFSS3.1 Vectrino Cross-Section 60mm above bed','FontSize',20)
cb = colorbar('northoutside');
caxis([0 250])
set(get(cb,'title'),'string','\bf\itMinutes After Low Tide','FontSize',20);

axes(ax(2))
m = length(FSS2.time)-6;
n = length(FSS2.bin15);
c = paruly(m);
t = FSS2.id:m;
for i = 1:n
plot(dist,log10(FSS2.bin15(i,:)),'Color',c(t(i),:),'linewidth',2,...
    'Marker','o','MarkerFaceColor',[1 1 1],'MarkerSize',6)
hold on
end
grid on
set(ax(2),'XTickLabel',[])
text(0.055,0.92,'\bf\itVP1','units','normalized','FontSize',20)
text(0.6,0.92,'\bf\itVP2','units','normalized','FontSize',20)
text(0.9,0.92,'\bf\itVP3','units','normalized','FontSize',20)
ylabel(['\bf\it' '$Log_{10}\left(\hat{\epsilon}_{TKE}\right)$'],'interpreter','latex','FontSize',18)
line([0 0],[-2 4],'Color',[0 0 0],'linewidth',2)

axes(ax(3))
m = length(FSS2.time)-11;
n = length(FSS1.bin15);
c = paruly(m);
t = FSS1.id:m;
for i = 1:n
plot(dist,log10(FSS1.bin15(i,:)),'Color',c(t(i),:),'linewidth',2,...
    'Marker','o','MarkerFaceColor',[1 1 1],'MarkerSize',6)
hold on
end
grid on
text(0.055,0.92,'\bf\itVP1','units','normalized','FontSize',20)
text(0.6,0.92,'\bf\itVP2','units','normalized','FontSize',20)
text(0.9,0.92,'\bf\itVP3','units','normalized','FontSize',20)
xlabel('\bf\itDistance Along Cross-Shore Transect (cm)','FontSize',20)
line([0 0],[-2 4],'Color',[0 0 0],'linewidth',2)

%global plot adjustments
set(ax,'YLim',[-1.5 4],'YTick',-1.5:1:4,...
    'XLim',[-13 23],'XTick',-13:3:23,'GridLineStyle',':',...
    'FontSize',20,'FontName','Arial')
set(ax,'LineWidth',1.5)
set(ax(1),'position',[0.125 0.63 0.79 0.25])
set(ax(2),'position',[0.125 0.355 0.79 0.25])
set(ax(3),'position',[0.125 0.08 0.79 0.25])
set(cb,'linewidth',1.5)
if savefigs
    fpath = savefigdir;fname = [dname '_compare'];
    export_fig([fpath fname],'-pdf','-nocrop','-m1','-r900')
    disp(['Figure ' fname '.pdf saved'])
end

