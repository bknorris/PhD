load('C:\Users\bkn5\Projects\Mekong_W2015\Data\Vector\FSS\V5109_070315.mat')
start = datenum(2015,03,14,09,20,00);
stop = datenum(2015,03,14,19,00,00);
ind = find(ADV.datetime >= start & ADV.datetime <= stop);
ind2 = find(ADV.Sensor.Datetime >= start & ADV.Sensor.Datetime <= stop);
adv.u = ADV.U(ind,:);
adv.v = ADV.V(ind,:);
adv.p = ADV.Pres(ind,:);
adv.dtime = ADV.datetime(ind,:);
adv.lat = ADV.Metadata.inst_lat;
adv.lon = ADV.Metadata.inst_lon;
adv.name = 'V5109';
disp(['Compass bearing for ' adv.name ': ' num2str(nanmean(ADV.Sensor.Heading(ind2))) ' degrees'])


load('C:\Users\bkn5\Projects\Mekong_W2015\Data\Aquadopp\FSS\AD5116_9March2015.mat')
ind = find(aqdp.datenum >= start & aqdp.datenum <= stop);
a.u = aqdp.u(ind,:);
a.v = aqdp.v(ind,:);
a.dtime = aqdp.datenum(ind,:);
a.p = aqdp.pressure(ind,:);
a.lat = aqdp.metadata.lat;
a.lon = aqdp.metadata.lon;
a.name = 'AD5116';
a.date = datenum(aqdp.datenum(ind(1),:));


%for the ADV:
fs = 32;
avt = fs*60*10;
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

lat = a.lat;
lon = a.lon;

%depth average
u = nanmean(a.u,2);
v = nanmean(a.v,2);
fs = 8;%hz
avt = fs*60*10;

%plot pressure record for reference    
idx = [1 avt:avt:length(a.p)];
c = jet(length(idx));
%time average & plot quivers
for j = 1:length(idx)-1
    t = idx(j):idx(j+1);
    U = nanmean(u(t));
    V = nanmean(v(t));
    hold on
    quiver(lon,lat,U,V,0.001,'Color',c(j,:))
    text(lon,lat,a.name)
end

% load('C:\Users\bkn5\Projects\Mekong_W2015\Data\Vector\F2F2\VC01_050315.mat')
% start = datenum(2015,03,05,09,00,00);
% stop = datenum(2015,03,05,18,00,00);
% ind = find(ADV.datetime >= start & ADV.datetime <= stop);
% ind2 = find(ADV.Sensor.Datetime >= start & ADV.Sensor.Datetime <= stop);
% adv.u = ADV.U(ind,:);
% adv.v = ADV.V(ind,:);
% adv.p = ADV.Pres(ind,:);
% adv.dtime = ADV.datetime(ind,:);
% adv.lat = ADV.Metadata.inst_lat;
% adv.lon = ADV.Metadata.inst_lon;
% adv.name = 'V1';
% disp(['Compass bearing for ' adv.name ': ' num2str(nanmean(ADV.Sensor.Heading(ind2))) ' degrees'])
