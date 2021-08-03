%Plot wind speed and direction during the HTA-VTA experiment
clear
close all
load('d:\Projects\Mekong_W2015\Data\Weather\Mekong2015_metstation.mat')
time = met.datetime;time2 = (time);
spd = met.windspd;
dir = met.winddir+180;
%average spd,dir
fs = 0.0167;
avt = 60*120;
win = round(fs*avt);
ind = [1 win:win:length(spd)];
n = length(ind);
speed = zeros(n-1,1);
dirc = zeros(n-1,1);
time3 = zeros(n-1,1);
for i = 1:n-1;
    speed(i) = mean(spd(ind(i):ind(i+1)));
    dirc(i) = mean(dir(ind(i):ind(i+1)));
    time3(i) = time2(round(median(ind(i):ind(i+1))));
end
dirc(dirc>360) = NaN;
tmax = max(time3);tmin = min(time3);
%times
step1 = datenum(0,0,0,20,0,0);
step2 = datenum(0,0,0,10,0,0);
%experiment gaps
t1 = datenum(2015,03,07,13,36,00);e1 = datenum(2015,03,07,17,10,00);
t2 = datenum(2015,03,08,14,15,00);e2 = datenum(2015,03,08,19,00,00);
t3 = datenum(2015,03,10,14,45,00);e3 = datenum(2015,03,10,16,40,00);
t4 = datenum(2015,03,14,4,40,00);e4 = datenum(2015,03,14,10,40,00);

%plot routine
f1 = figure;
set(f1,'PaperOrientation','portrait',...
    'position',[500 100   1000   300]);
hold on
p = patch([t1 e1 e1 t1],[0 0 12 12],[.9 .9 .9]);set(p,'EdgeColor','none')
p = patch([t2 e2 e2 t2],[0 0 12 12],[.9 .9 .9]);set(p,'EdgeColor','none')
p = patch([t3 e3 e3 t3],[0 0 12 12],[.9 .9 .9]);set(p,'EdgeColor','none')
p = patch([t4 e4 e4 t4],[0 0 12 12],[.9 .9 .9]);set(p,'EdgeColor','none')

b(1) = line(time3,speed,...
    'linewidth',1.5,...
    'Color',[0.4 0.4 0.4]);
%text labels
text(68,10,'HTA')
text(73,10,'VTA')
set(b,'linestyle','--')
ax1 = gca;
set(ax1,'XColor','k');
set(ax1,'YColor',[0.4 0.4 0.4]);
ax1_pos = get(ax1,'position'); % position of first axes
ax2 = axes('Position',ax1_pos,...
    'XAxisLocation','top',...
    'YAxisLocation','right',...
    'Color','none');
line(time3,dirc,...
    'Parent',ax2,...
    'Color','k',...
    'linewidth',1.5)
set(ax2,'xticklabel',[],...
    'ytick',0:90:360,...
    'ylim',[0 360],...
    'xlim',[tmin tmax])
set(ax1,'ylim',[0 12],'xlim',[tmin tmax])
datetick(ax1,'x','dd','keepticks','keeplimits')
xlabel(ax1,'Day in March, 2015')
ylabel(ax1,'Wind Speed (ms^-^1)')
ylabel(ax2,'Wind Direction (deg)')
prettyfigures('text',13,'labels',14,'box',0)
figdir = 'f:\GradSchool\DataAnalysis\Paper2\WorkingFigures\Environmental\';
export_fig([figdir 'WeatherTimeseries'],'-pdf','-nocrop')

