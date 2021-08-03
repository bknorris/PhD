%Testbed for wavelet analysis with velocity and bottom trace

clear
close all
% savefigdir = 'c:\Users\Bnorr\Documents\GradSchool\DataAnalysis\Paper3\WorkingFigures\';
load('d:\Mekong_W2015\DataAnalysis\Paper3\BottomTrack\08-03-15\HTA2_bdtrace.mat')
bd = vpro1.bdist;bdt = vpro1.time;clear vpro1 vpro2 vpro3
load('d:\Mekong_W2015\DataAnalysis\Paper3\VPs\8March2015_Vels.mat')
x = vpro1.x;y = vpro1.y;vt = vpro1.time;
clear vpro1 vpro2 vpro3
t1 = bdt(1);t2 = bdt(end);
id = find(vt>= t1 & vt<= t2);
x = nanmean(x(id,5:25),2);y = nanmean(y(id,5:25),2);vt = vt(id);
%rotate
rot = -22*pi/180;
x = x.*(ones(size(x))*cos(rot)) + ...
    y.*(ones(size(y))*sin(rot));
y = -y.*(ones(size(y))*sin(rot)) + ...
    x.*(ones(size(x))*cos(rot));
%interp to 10Hz
[vt,idx]=unique(vt);
x10 = interp1(vt,x(idx),bdt);
y10 = interp1(vt,y(idx),bdt);
bdh = 0.24-bd;

% bds = fastsmooth(bdh,floor(length(bdh)/20),1,1);
% brsd = bsxfun(@minus,(bdh-bds),mean(bdh-bds)); %mean removal & detrend
% xs = detrend(x10)t = bdt;

figure
p(1) = subplot(211);
plot(bdt,gradient(bdh),'k','linewidth',1.5),hold on
ylabel('BLE [m]')
p(2) = subplot(212);
plot(bdt,x10,'k','linewidth',1.5)
set(p(1),'xticklabel','')
ylabel('X-shore Velocity [m/s]')
xlabel('Time on March 8th, 2015')
datetickzoom('x','HH:MM:SS','keepticks','keeplimits')
linkaxes(p,'x')
clear p;
prettyfigures('text',13,'labels',14,'box',1)

x10 = detrend(x10);
bdh = detrend(bdh);
%Use a sample of the data
b = bdh(20*60*10:40*60*10);
e = x10(20*60*10:40*60*10);
ts = bdt(20*60*10:40*60*10);

dt = 1/10;
time = (dt:dt:length(ts)*dt)./60; %time in minutes for x-axis plotting
xt = [time' e];
yt = [time' gradient(b)];
[Rsq,period,scale,coi,sig95,Wxy,t,dt] = wtc(xt,yt,'Dj',(1/10),'S0',(1/60));

threshold = 0.9;
[m,n] = size(sig95);
events = zeros(m,n);
for k = 1:m
    events(k,:) = sig95(k,:) >= threshold;
end
%filter by coi
zid = zeros(m,n);
for k = 1:n
    zid(:,k) = period <= coi(k);
end
events1 = events.*zid;
eventGroups1 = bwlabel(events,8);

dt = 1/10;
time = (dt:dt:length(ts)*dt)./60; %time in minutes for x-axis plotting
xt = [time' e];
yt = [time' gradient(b)];
[Rsq,period,scale,coi,sig95,Wxy,t,dt] = wtc(xt,yt,'S0',10/60);

threshold = 0.9;
[m,n] = size(sig95);
events = zeros(m,n);
for k = 1:m
    events(k,:) = sig95(k,:) >= threshold;
end
%filter by coi
zid = zeros(m,n);
for k = 1:n
    zid(:,k) = period <= coi(k);
end
events = events.*zid;
eventGroups2 = bwlabel(events,8);
figure
subplot(121)
imagesc(events1)
subplot(122)
imagesc(events2)
figure
imagesc(eventGroups1-eventGroups2)
