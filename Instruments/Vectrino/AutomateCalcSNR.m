%Two routines: the first plots raw SNR and velocities from a given
%deployment as an example. The second loops through the VP data files and
%saves out the SNR as a separate file for future QC with the TKE
%dissipation rates. 
clear
maindir = 'D:\Projects\Mekong_F2014\Data\Vectrino\';
veldir = 'd:\Projects\Mekong_F2014\DataAnalysis\Paper2\';
folder = 'DPS';
file = 'VP1_021014.mat';
gmt2ict = datenum(0,0,0,6,0,0); %6hrs for 2014, 7hrs for 2015
load([maindir folder '\' file])

SNR = (VPRO.Data.Profiles_SNRBeam1+VPRO.Data.Profiles_SNRBeam3)./2;
rb = VPRO.Data.Profiles_Range;
time = VPRO.Data.Time;
step = (time(end)-time(1))./4;

clear VPRO
load([veldir 'DPSVels.mat'])
savefigdir = 'd:\Projects\Mekong_F2014\Figures\Paper2\';

%%%Quickplot SNR and velocity data
f1 = figure(1);
set(f1,'PaperOrientation','portrait',...
        'position',[400 200   1000   500]);
    set(gcf,'color','w','PaperPositionMode','auto')
p(1) = subplot(121);
imagesc(time,rb,SNR')
cb = colorbar;
set(gca,'Xlim',[time(1) time(end)],'XTick',time(1):step:time(end),...
    'YLim',[rb(1) rb(end)],'YDir','normal')
caxis([15 60])
datetick('x','HH:MM','keepticks','keeplimits')
ylabel(cb,'Signal to Noise Ratio (dB)')
xlabel(['Time on ' datestr(time(1),'dd/mm/yy')])
p(2) = subplot(122);
imagesc(time,rb,x')
cb = colorbar;
set(gca,'Xlim',[time(1) time(end)],'XTick',time(1):step:time(end),...
    'YLim',[rb(1) rb(end)],'YDir','normal')
caxis([-0.02 0.02])
datetick('x','HH:MM','keepticks','keeplimits')
ylabel(cb,'Velocity (m/s)')
xlabel(['Time on ' datestr(time(1),'dd/mm/yy')])

