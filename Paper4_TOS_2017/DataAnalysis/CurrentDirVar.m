load('D:\Projects\Mekong_W2015\DataAnalysis\Paper1\HTA_1Vels.mat')
start = datenum(2015,03,07,14,24,00); stop = datenum(2015,03,07,15,21,00);
id = find(dat.vpro2.time >= start & dat.vpro2.time <= stop);
time = dat.vpro2.time(id,:);
x = mean(dat.vpro2.x(id,1:5),2);
y = mean(dat.vpro2.y(id,1:5),2);
[th,r] = cart2compass(y,x);
% ids = find(th > 0 & th < 180);
% th(ids) = th(ids)+360;
avt = 50*60*10;
ind = [1 avt:avt:length(x)];
figure
% for i = 1:length(ind)-1
% v = var(th(ind(i):ind(i+1)));
% time2 = time(ind(i));
% plot(time2,v,'*r'),hold on
% end
plot(time,th,'*r')
datetickzoom('x','HH:MM','keepticks','keeplimits')
xlabel(['Time on ' datestr(time(1),'dd-mm-yyyy')])
ylabel('Dir (degrees)')
title('Variance in direction (10 min), HTA1 VP2 bins 1:5')
ylabel('Var(dir)')