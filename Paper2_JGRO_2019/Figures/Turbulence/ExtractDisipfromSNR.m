%Quick plot SNR and extracted TKE values

clear
%%%% EXAMPLE PLOTS %%%%
load('D:\Projects\Mekong_W2015\DataAnalysis\Paper2\HTA_4Vels.mat')
tt = datenum(2015,03,10,16,00,00);
tmp = abs(dat.vpro2.time-tt);
[~,idx] = min(tmp);
snr = dat.vpro2.SNR1(idx,:);
rb = dat.vpro2.rb;
savefigdir = 'd:\Projects\Mekong_W2015\Figures\Paper2\';

f1 = figure(1);
set(f1,'PaperOrientation','portrait',...
    'position',[400 200   600   500]);
set(gcf,'color','w','PaperPositionMode','auto')
plot(snr,rb,'LineWidth',1.5,'Color','k')
set(gca,'ydir','reverse')
ylabel('Range (m)')
xlabel('SNR (dB)')
title('HTA Day 1 VP2 SNR at 15:50:00')
export_fig([savefigdir 'SNR_profile'],'-jpeg','-nocrop')

time = dat.vpro2.time;snr = dat.vpro2.SNR1;
f2 = figure(2);
set(f2,'PaperOrientation','portrait',...
    'position',[400 200   600   500]);
set(gcf,'color','w','PaperPositionMode','auto')
imagesc(time,rb,snr'),datetick('x','HH:MM','keepticks','keeplimits')
caxis([40 60])
cb = colorbar;ylabel(cb,'dB')
set(gca,'ydir','reverse')
ylabel('Range (m)')
xlabel('Time on 07/03/2015')
title('HTA Day 1 VP2 SNR')
export_fig([savefigdir 'SNR_colorts'],'-jpeg','-nocrop')

vel = dat.vpro2.beam1;
f3 = figure(3);
set(f3,'PaperOrientation','portrait',...
    'position',[400 200   600   500]);
set(gcf,'color','w','PaperPositionMode','auto')
imagesc(time,rb,vel'),datetick('x','HH:MM','keepticks','keeplimits')
caxis([-0.1 0.1])
cb = colorbar;ylabel(cb,'m/s')
set(gca,'ydir','reverse')
ylabel('Range (m)')
xlabel('Time on 07/03/2015')
title('HTA Day 1 VP2 beam1 velocities')
export_fig([savefigdir 'Vel_colorts'],'-jpeg','-nocrop')
clear dat

%%%% Main Routine %%%%
%Load TKE file and corresponding SNR file
tkedir = 'd:\Projects\Mekong_W2015\DataAnalysis\Paper2\TKE\';
snrdir = 'd:\Projects\Mekong_W2015\DataAnalysis\Paper2\SNR\';
load([tkedir 'HTA_4TKE.mat'])
load([snrdir 'HTA_4SNRavgs.mat'])

epsnr = zeros(length(Stat.vpro2.beam1.E),1);
%average SNR together
signal = (SNR.vpro2.beam1+SNR.vpro2.beam2+SNR.vpro2.beam3+SNR.vpro2.beam4)./4;
eps = (Stat.vpro2.beam1.E+Stat.vpro2.beam2.E+Stat.vpro2.beam3.E+Stat.vpro2.beam4.E)./4;
for i = 8:length(signal)
    ok = 0;
    while ok < 1
        [~,id] = max(signal(i,:));
        if isnan(eps(id,i)) || eps(id,i) == 0
            signal(i,id) = 0;
        else
            ok = 1;
            epsnr(i) = eps(id,i);
        end
    end
end
epsold = eps(7,:);
time = SNR.vpro2.time;time(time < 1E5) = NaN;

f4 = figure(4);
set(f4,'PaperOrientation','portrait',...
    'position',[400 200   600   500]);
set(gcf,'color','w','PaperPositionMode','auto')
p(1) = plot(time,epsold,'LineWidth',1.5,'Color','k'); hold on
p(2) = plot(time,epsnr,'LineWidth',1.5,'Color','r');
leg = legend(p,{'Bin 15 Method';'Max SNR Method'});set(leg,'box','off')
set(gca,'Xlim',[time(1) time(end)])
datetick('x','HH:MM','keepticks','keeplimits')
ylabel('Range (m)')
xlabel('Time on 07/03/2015')
title('HTA Day 1 VP2 Dissipaton Rate Method Compare')
export_fig([savefigdir 'TKEmethodcmp'],'-jpeg','-nocrop')

f5 = figure(5);
set(f5,'PaperOrientation','portrait',...
    'position',[400 200   600   500]);
set(gcf,'color','w','PaperPositionMode','auto')
imagesc(time,rb,Stat.vpro2.beam1.E)
set(gca,'Xlim',[time(1) time(end)])
datetick('x','HH:MM','keepticks','keeplimits')
caxis([0.0005 0.005])
cb = colorbar;ylabel(cb,'m^2s^-^3')
set(gca,'ydir','reverse')
ylabel('Range (m)')
xlabel('Time on 07/03/2015')
title('HTA Day 1 VP2 Averaged Dissipation Rate')
export_fig([savefigdir 'TKEcolorts'],'-jpeg','-nocrop')
