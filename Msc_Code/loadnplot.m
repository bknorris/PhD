clear

%set these first
start = datenum(2014,06,18,14,53,00);
finish = datenum(2014,06,18,15,26,00);

%% Import SLOBS
fname = 'red_snd_10_100_1000_SLOBS_dat.csv';
fid = fopen(fname);
str = fgetl(fid);
header = textscan(str,'%s','Delimiter',',');
fields = header{1};
FRMT = '%s%n%n%n';
thedata = textscan(fid,FRMT,'delimiter',',');
slobs = cell2struct(thedata,fields,2);

%Convert time into serial date no
timestr = thedata{1};
thedate = {'18-Jun-2014 '}; %date of the experimental run
datefrmt = 'dd-mmm-yyyy HH:MM:SS';
datestr1 = repmat(thedate,length(timestr),1);
datestr2 = strcat(datestr1,timestr);
datetime = datenum(datestr2,datefrmt);
slobs = rmfield(slobs,'Timecode');
slobs.datetime = datetime;

%Crop data to start/stop times
ind = find(datetime >= start & datetime <= finish);
slobs.exptime = datetime(ind);
% ind2 = find(slobs.Gain < 20);slobs.FTU(ind2) = NaN; %find obscured values
slobs.expftu = slobs.FTU(ind);

%Average Data
int = 5; %sample interval (in seconds)
aveint = 60; %average interval (in seconds)
nsamp = numel(slobs.expftu);
subsamp = aveint/int;
i = 1:subsamp:nsamp;
for ii = 1:length(i)-1
    slobs.average(ii) = sum(slobs.expftu(i(ii):i(ii+1)))./subsamp;
    slobs.stdev(ii) = std(slobs.expftu(i(ii):i(ii+1)));
end

%% Import PAR Sensor
fname = 'red_snd_10_100_1000_PAR_dat.txt';
fid = fopen(fname);
str = fgetl(fid);
header = textscan(str,'%s','Delimiter',',');
fields = {'dummy','datetime','val'};
FRMT = '%n%s%n';
thedata = textscan(fid,FRMT,'delimiter',',');
par = cell2struct(thedata,fields,2);
par = rmfield(par,'dummy');

%Convert time into serial date no
datetime = cellfun(@datenum,par.datetime);

%Crop data to start/stop times
ind = find(datetime >= start & datetime <= finish);
par.exptime = datetime(ind);
expval = par.val(ind);
par.expval = (expval/4.6); %convert to W/m^2

%Average Data
int = 1; %sample interval (in seconds)
aveint = 60; %average interval (in seconds)
nsamp = numel(par.expval);
subsamp = aveint/int;
i = 1:subsamp:nsamp;
for ii = 1:length(i)-1
    par.average(ii) = sum(par.expval(i(ii):i(ii+1)))./subsamp;
    par.stdev(ii) = std(par.expval(i(ii):i(ii+1)));
end
par.average(33) = (par.average(31)+par.average(32))/2.5;
par.stdev(33) = 0;

%% Plots
depth = 0:5:160; %step interval (in cm)

m1 = figure('PaperOrientation','portrait',...
    'position',[400   100   800   800]);

axes('position',[0.1 0.1 0.6 0.8])
t1 = line(slobs.average,depth);
hold on
h1 = herrorbar(slobs.average,depth,slobs.stdev);
hold on
set(t1,'LineWidth',1,'Color','b')
haxes1 = gca;
haxes1_pos = get(haxes1,'Position'); % store position of first axes
haxes2 = axes('Position',haxes1_pos,...
              'XAxisLocation','top',...
              'YAxisLocation','right',...
              'Color','none');
t2 = line(par.average,depth,'Parent',haxes2);
hold on
h2 = herrorbar(par.average,depth,par.stdev);
hold on
set(t2,'LineWidth',1,'Color','r')
set(h2,'LineWidth',1,'Color','r')
set(haxes1,...
    'YDir','Reverse',...
    'YGrid','on',...
    'TickDir','out',...
    'XAxisLocation','bottom',...
    'XColor',[.2 .2 .2], ...
    'YColor',[.2 .2 .2], ...
    'LineWidth',1,...
    'TickLength',[.02 .02],...
    'YMinorTick','on',...
    'YTick',[0:10:160],...
    'YLim', [0 160])
set(haxes2,... 
    'YDir','Reverse',...
    'YTickLabel',[],...
    'LineWidth',1,...
    'TickDir','out',...
    'XColor',[.2 .2 .2], ...
    'YColor',[.2 .2 .2], ...    
    'TickLength',[.01 .01],...
    'YMinorTick','on',...
    'YTick',[0:10:160],...
    'YLim', [0 160])
leg = legend([t1 t2],{'\bf\itSLOBS' '\bf\itPAR Sensor'},'location','SouthEast');
ylabel(haxes1,'\itDepth (cm)')
xlabel(haxes1,'\itFTU')
xlabel(haxes2,'\itW/m^2')
title('\bfPAR Attenuation and Turbidity Variance by Depth, 1000mg/L Assumed Concentration')

ylim=get(gca,'YLim');
xlim=get(gca,'XLim');
text(xlim(2)+2,depth(1),'\bf1 Minute Average (FTU)',...
   'VerticalAlignment','bottom',...
   'HorizontalAlignment','left')
for i = 2:length(slobs.average)
    text(xlim(2)+5,depth(i),num2str(slobs.average(i)),...
   'VerticalAlignment','bottom',...
   'HorizontalAlignment','left')
end
box on
set(gcf, 'PaperPositionMode', 'auto');

% print '-dpng' 'red_snd_1000mg.png'

save D:\Projects\PARSensor\Sand\red_snd_1000_PAR par
save D:\Projects\PARSensor\Sand\red_snd_1000_SLOBS slobs

pause(5)
close all

