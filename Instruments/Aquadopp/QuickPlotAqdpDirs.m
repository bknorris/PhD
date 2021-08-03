%check the aquadopp directions with the vector from F2F2 and F2F3
clear
aqdpdir = 'D:\Projects\Mekong_F2014\Data\Aquadopp\DPS\';
vecdir = 'D:\Projects\Mekong_F2014\Data\Vector\DPS\';

aqdps = dir([aqdpdir '*.mat']);
%sort dir by date, use only the most recent files
dd = strfind({aqdps.date},'29-Aug-2016'); %date the files were reprocessed
id = find(not(cellfun('isempty', dd)));
name = {'AD5116';'AD5117';'LR4'};
%run through all the aquadopps
start = datenum(2014,10,02,11,00,00);
stop = datenum(2014,10,03,02,00,00);
count = 1;
for i = 1:length(id)
    iid = ['a' num2str(count)];
    disp(['Loading ' aqdps(id(i)).name])
    load([aqdpdir aqdps(id(i)).name])
    ind = find(aqdp.datenum >= start & aqdp.datenum <= stop);
    A.(iid).u = aqdp.u(ind,:);
    A.(iid).v = aqdp.v(ind,:);
    A.(iid).dtime = aqdp.datenum(ind,:);
    A.(iid).p = aqdp.pressure(ind,:);
    A.(iid).lat = aqdp.metadata.lat;
    A.(iid).lon = aqdp.metadata.lon;
    A.(iid).name = name(count);
    A.(iid).date = datenum(aqdp.datenum(ind(1),:));
    disp(['Compass bearing for ' A.(iid).name{:} ': ' num2str(nanmean(aqdp.heading(ind))) ' degrees'])
    clear aqdp
    count = count+1;
end

%load the vector
vec = dir([vecdir '*.mat']);
disp(['Loading ' vec(1).name])
load([vecdir vec(1).name])
ind = find(ADV.datetime >= start & ADV.datetime <= stop);
ind2 = find(ADV.Sensor.Datetime >= start & ADV.Sensor.Datetime <= stop);
adv.u = ADV.U(ind,:);
adv.v = ADV.V(ind,:);
adv.p = ADV.Pres(ind,:);
adv.dtime = ADV.datetime(ind,:);
adv.lat = ADV.Metadata.inst_lat;
adv.lon = ADV.Metadata.inst_lon;
adv.name = 'VC101';
disp(['Compass bearing for ' adv.name ': ' num2str(nanmean(ADV.Sensor.Heading(ind2))) ' degrees'])
clear ADV

% %load met data
% 
% metdir = 'c:\Users\bkn5\Projects\Mekong_W2015\Data\Weather\';
% load([metdir 'Mekong2015_metstation.mat'])
% ind = find(met.datetime >= start & met.datetime <= stop);
% m.t = met.datetime(ind,:);
% wdir = abs(met.winddir(ind,:)+180); %to make this TO not FROM
% wspd = met.windspd(ind,:);
% [m.u,m.v] = cmgspd2uv(wspd,wdir);
% m.lat = met.metadata.latitude;
% m.lon = met.metadata.longitude;
% m.name = 'MET';

fn = fieldnames(A);
f1 = figure(1);
set(f1,'PaperOrientation','landscape',...
    'position',[200 200   1200  600]);
set(gcf,'color','w','PaperPositionMode','auto')    

subplot(121)
%for the aquadopps
intv = 10; %minutes
fs = 8;%hz
avt = fs*60*intv;

%plot pressure record for reference    
idx = [1 avt:avt:length(A.a1.p)];
c = jet(length(idx));
for j = 1:length(idx)-1
    t = idx(j):idx(j+1);
    plot(A.a1.dtime(t),A.a1.p(t),'Color',c(j,:)), hold on
end
title(datestr(A.a1.date,'dd-mm-yyyy'))
datetick('x','HH:MM:SS','keepticks','keeplimits')
hold off
 
subplot(122)
for ii = 1:3 %cycle through instruments
    iid = fn{ii};
    if ii == 3
        fs = 1;
        avt = fs*60*intv;
        idx = [1 avt:avt:length(A.a3.p)];
    end
    %need to depth then time average, then calculate quivers at a given
    %time interval... say 10 minutes.
    lat = A.(iid).lat;
    lon = A.(iid).lon;
    
    %depth average
    u = nanmean(A.(iid).u,2);
    v = nanmean(A.(iid).v,2);

    %time average & plot quivers
    for j = 1:length(idx)-1
        t = idx(j):idx(j+1);
        U = nanmean(u(t));
        V = nanmean(v(t));
        hold on
        quiver(lon,lat,U,V,0.001,'Color',c(j,:))
        text(lon,lat,A.(iid).name)
    end
end

%for the ADV:
fs = 32;
avt = fs*60*intv;
idx = [1 avt:avt:length(adv.u)];
c = jet(length(idx));
figure(2)
for j = 1:length(idx)-1
    t = idx(j):idx(j+1);
    plot(adv.dtime(t),adv.p(t),'Color',c(j,:)), hold on
end
title(['ADV Time ' datestr(adv.dtime(1),'dd-mm-yyyy')])
datetick('x','HH:MM:SS','keepticks','keeplimits')
hold off
figure(1)
for j = 1:length(idx)-1
    t = idx(j):idx(j+1);
    advU = nanmean(adv.u(t));
    advV = nanmean(adv.v(t));
    quiver(adv.lon,adv.lat,advU,advV,0.001,'Color',c(j,:)), hold on
    text(adv.lon,adv.lat,adv.name)
end

% %for the met station
% fs = 0.1;
% avt = fs*60*intv;
% idx = [1 avt:avt:length(m.u)];
% c = jet(length(idx));
% for j = 1:length(idx)-1
%     t = idx(j):idx(j+1);
%     mu = nanmean(m.u(t));
%     mv = nanmean(m.v(t));
%     quiver(m.lon,m.lat,mu,mv,0.0001,'Color',c(j,:)), hold on
%     text(m.lon,m.lat,m.name)
% end
% hold off
% savedir = 'c:\users\bkn5\desktop\rotations\';
% set(gca,'Xlim',[106.2436 106.245],'YLim',[9.4905 9.4935],'box','on')
% title(datestr(A.a1.date,'dd-mm-yyyy'))
% ylabel('Latitude'),xlabel('Longitude')
% colormap(jet)
% export_fig([savedir 'Dep' datestr(A.a1.date,'dd-mm-yyyy')],'-png')
% 
