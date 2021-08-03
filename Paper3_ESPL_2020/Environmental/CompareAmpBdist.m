%Plot bottom surface track, backscatter and compare.
clear
close all

%F2F and F2F2 both have instruments deployed at 80mm, close enough to see
%the bed in the Backscatter.
dir1 = 'e:\Mekong_W2015\DataAnalysis\Paper3\VPs\06-03-15\';
files = {'6March2015_Vels.mat';'6March2015_Sen.mat'};
%VP2
vars = whos('-file',[dir1 files{1}]);vars = {vars.name};
V = matfile([dir1 files{1}]);
S = matfile([dir1 files{2}]);

i = 1;
vels = V.(vars{i});
sen = S.(vars{i});
vph = 0.063;

start = datenum('06-Mar-2015 13:37:40');
stop = datenum('06-Mar-2015 16:18:58');
rb = vels.rb;h = vph-rb;
bdtime = vels.bdtime;
id = find(bdtime >= start & bdtime <= stop);
bdtime = bdtime(id);bdist = vph-vels.bdist(id);
time = vels.time;

%extract bottom distance
bd = my_running_median(bdist,512);

%compare to backscatter
id = find(time >= start & time <= stop);
time = time(id);
amp = sen.Amp3(id,:);
bd50 = spline(bdtime,bd,time);
%find max(amp)
nsamp = length(time);
amax = zeros(nsamp,1);
for ii = 1:nsamp
    id = find(h<=bd50(ii)+0.002 & h>=bd50(ii)-0.002);
    if isempty(id)
        continue
    end
    [~,aid] = max(amp(ii,id));
    amax(ii) = h(id(aid));
end
%plot
f1 = figure(1);
set(f1,'PaperOrientation','portrait',...
    'position',[400 100   500  350],...
    'renderer','painters');
nudge = datenum(0,0,0,0,10,0);
sp(1) = subplot(211);
plot(downsample(bdtime,100),downsample(bdist,100),'k'),hold on
plot(bdtime,bd,'r',...
    'linewidth',1.5)
sp(2) = subplot(212);
cmap = brewermap(50,'YlOrRd');
colormap(cmap);
imagesc(time,h,amp'),hold on
plot(time,fastsmooth(amax,50,3,1),'w')
plot(time,bd50,'k','linewidth',2)
%plot adjustments
set(sp(1),'xlim',[bdtime(1)+nudge bdtime(end)-nudge],...
    'ylim',[-0.01 0.05],...
    'xticklabel',[],...
    'position',[0.15 0.58 0.65 0.35])
set(sp(2),'xlim',[bdtime(1)+nudge bdtime(end)-nudge],...
    'ylim',[-0.01 0.025],...
    'ydir','normal',...
    'position',[0.15 0.14 0.65 0.35])
caxis([-20 -5])
datetick('x','HH:MM','keepticks','keeplimits')
xlabel(sp(2),['Time on ' datestr(time(1),'dd-mm-yy')])
suplabel('Height above bottom [m]','y')
cb = colorbar;
set(cb,'position',[0.82 0.14 0.04 0.35],...
    'ytick',-20:5:-5,...
    'linewidth',1.5,...
    'tickdir','out')
ylabel(cb,'Amplitude (dB)')
prettyfigures('text',12,'labels',13,'box',1)
sfdir = 'd:\GradSchool\DataAnalysis\Paper3\Figures\';
export_fig([sfdir 'BdAmpCompare'],'-pdf')