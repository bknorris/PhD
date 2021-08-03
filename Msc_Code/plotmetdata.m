%plot met data
clear
load metdata.mat

startt = datenum(2015,03,07,08,00,00);
endt = datenum(2015,03,11,09,00,00);
fname = 'FineScaleStudy3';
ind = find(met.datetime >= startt & met.datetime <= endt);

%met variables to plot:
%Wind speed, atmosphere pres, wind direction, cross/along shore wind speed, temperature
%first crop variables:

time = met.datetime(ind);
pres = met.pressure(ind);
wndspd = met.windspd(ind);
temp = met.temp(ind);
winddir = met.winddir(ind);

[windu,windv] = cmgspd2uv(wndspd,winddir);

f1 = figure;
set(f1,'PaperOrientation','portrait',...
    'position',[400 200   1400   700]);
ax = tight_subplot(5,1,0.1,[0.1 0.05],[0.1 0.05]);
axes(ax(1))
plot(time,wndspd,'b')
ylabel('m s-1')
leg = legend('Wind Speed');
set(leg,'position',[0.87 0.95 0.1 0.01])
set(gca,'Ylim',[0 15])
grid on

axes(ax(2))
plot(time,winddir,'+b')
ylabel('degrees')
leg = legend('Wind Dir');
set(leg,'position',[0.87 0.75 0.1 0.01])
set(gca,'Ylim',[0 360],'YTick',0:90:360)
grid on

axes(ax(3))
plot(time,pres,'-k')
ylabel('mbar')
leg = legend('Atmos P');
set(leg,'position',[0.87 0.57 0.1 0.01])
set(gca,'Ylim',[1005 1020])
grid on

axes(ax(4))
plot(time,windu,'-b')
hold on
plot(time,windv,'-r')
ylabel('m s-1')
leg = legend('Along Shore','Cross Shore');
set(leg,'position',[0.87 0.37 0.1 0.05])
set(gca,'Ylim',[-10 10])
grid on

axes(ax(5))
plot(time,temp,'r')
ylabel('C')
leg = legend('Temp');
set(leg,'position',[0.87 0.19 0.1 0.01])
set(gca,'Ylim',[26 30])
datetickzoom('x','dd HH:MM:SS','keepticks','keeplimits')
xlabel('dd HH:MM:SS')
grid on
set(ax(5),'Xlimmode','auto','GridLineStyle',':','ticklength',[0.001 0.01])
set([ax(1) ax(2) ax(3) ax(4)],'XTickLabel',[],'Xlimmode','auto',...
            'GridLineStyle',':','ticklength',[0.001 0.01])
set(gcf,'color','w','PaperPositionMode','auto')
fpath = 'C:\Users\bkn5\Projects\Mekong_W2015\Figures\Met\';
export_fig([fpath fname],'-png','-m1','-r900','-opengl')
        