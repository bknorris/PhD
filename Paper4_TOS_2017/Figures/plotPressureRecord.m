%Plot the overview figures from the HTA and VTA. This is the pressure
%signal from these deployments. 

figdir = 'C:\Users\bkn5\Projects\Mekong_W2015\Figures\Environment\';
vdir = 'C:\Users\bkn5\Projects\Mekong_W2015\Data\Vector\SW_Mudflat\';

%the HTA took place over three days; the mudflat ADV is in two separate
%files over the course of that experiment. First, combine the pressure
%records to plot the pressure record for the whole experiment.
load([vdir 'V5108_030315.mat'])
ft1 = datenum(2015,03,07,00,00,00);fe1 = ADV.datetime(end);
ind = find(ADV.datetime >= ft1 & ADV.datetime <= fe1);
advtime1 = ADV.datetime(ind);advpres1 = smooth(ADV.Pres(ind),1000,'moving');
advdepth1 = advpres1+(ADV.Metadata.pressure_sensor_height/1000);
clear ADV
load([vdir 'V5108_080315.mat'])
ft2 = ADV.datetime(1);fe2 = datenum(2015,03,10,20,30,00);
ind = find(ADV.datetime >= ft2 & ADV.datetime <= fe2);
advtime2 = ADV.datetime(ind);advpres2 = smooth(ADV.Pres(ind),1000,'moving');
advdepth2 = advpres2+(ADV.Metadata.pressure_sensor_height/1000);
clear ADV

%combine
advtime = [advtime1; advtime2];advdepth = [advpres1; advpres2];

%experimental times, save into structure HTA:
t1 = datenum(2015,03,07,13,36,00);e1 = datenum(2015,03,07,17,10,00);HTA.times.t1 = t1;HTA.times.e1 = e1;
t2 = datenum(2015,03,08,14,15,00);e2 = datenum(2015,03,08,19,00,00);HTA.times.t2 = t2;HTA.times.e2 = e2;
t3 = datenum(2015,03,10,14,45,00);e3 = datenum(2015,03,10,16,30,00);HTA.times.t3 = t3;HTA.times.e3 = e3;

%plot data range figure
f1 = figure(1);
set(f1,'PaperOrientation','portrait',...
    'position',[400 100   1000 600]);
set(gcf,'color','w','PaperPositionMode','auto')
hold on
p = patch([t1 e1 e1 t1],[0 0 2 2],[.9 .9 .9]);set(p,'EdgeColor','none')
p = patch([t2 e2 e2 t2],[0 0 2 2],[.9 .9 .9]);set(p,'EdgeColor','none')
p = patch([t3 e3 e3 t3],[0 0 2 2],[.9 .9 .9]);set(p,'EdgeColor','none')
plot(advtime,advdepth,'k','LineWidth',1.5)
hold off
set(gca,'Ylim',[0 1.4],'Xlim',[ft1 fe2],'box','on','LineWidth',1.5,'FontSize',12)
datetick('x','mm/dd HH:MM','keepticks','keeplimits')
xlabel('\bf\itDays in 03/2015')
ylabel('\bf\itDepth Above Bed (m)'), hold off
clear ADV

%save figs
fpath = figdir;fname = 'HTApresRecord';
export_fig([fpath fname],'-pdf')
disp(['Figure ' fname '.pdf saved'])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Now, do the same for the mudflat Vector from the VTA
vdir = 'C:\Users\bkn5\Projects\Mekong_W2015\Data\Vector\DPS2\';
t1 = datenum(2015,03,14,04,10,00);e1 = datenum(2015,03,14,15,00,00);
t2 = datenum(2015,03,14,07,00,00);e2 = datenum(2015,03,14,11,00,00);
%t2 and e2 are the times when the top Vectrino of the VTA is underwater
VTA.times.t1 = t1;VTA.times.e1 = e1;VTA.times.t2 = t2;VTA.times.e2 = e2;

load([vdir 'V5109_130315.mat'])
ind = find(ADV.datetime >= t1 & ADV.datetime <= e1);
advtime = ADV.datetime(ind);advpres = smooth(ADV.Pres(ind),1000,'moving');
advdepth = advpres+(ADV.Metadata.pressure_sensor_height/1000);

%plot data range figure
f2 = figure(2);
set(f2,'PaperOrientation','portrait',...
    'position',[400 100   800 600]);
set(gcf,'color','w','PaperPositionMode','auto')
hold on
c = [0.7 0.7 0.7;0.5 0.5 0.5;0 0 0];
p = patch([t2 e2 e2 t2],[0 0 2 2],[.95 .95 .95]);set(p,'EdgeColor','none')
plot(advtime,advdepth,'k','LineWidth',1.5)
hold off
set(gca,'Ylim',[0 1.8],'Xlim',[t1 e1],'box','on','LineWidth',1.5,'FontSize',12)
datetick('x','HH:MM:SS','keepticks','keeplimits')
xlabel('\bf\itTime on 14/03/2015')
ylabel('\bf\itDepth Above Bed (m)'), hold off
clear ADV

%save figs
fpath = figdir;fname = 'VTApresRecord';
export_fig([fpath fname],'-pdf')
disp(['Figure ' fname '.pdf saved'])

clearvars -except VTA HTA
