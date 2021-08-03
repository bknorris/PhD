%plot the variance of depth and raw depth profiles for selected times from
%the various HTA codes (HTAday1, HTAday2, HTAday3). 

clear
load('C:\Users\bkn5\Projects\Mekong_W2015\Data\Aquadopp\FSS\AD5116_9March2015.mat')
%might need to do some QC on the data before using it... low correlations
%in the upper bins. 

%time average vertical velocity:
start = datenum(2015,03,07,15,30,00);stop = datenum(2015,03,07,15,40,00);
heading = 110; %cross shore angle
idx = find(aqdp.datenum >= start & aqdp.datenum <= stop);
h = 63.2; %inst hab (cm)
ph = 43.5; %pneumatophore height (cm)
badbins = 1:4; %some bins contain many low correlations. Designate them.
nCells = aqdp.nCells;
b1 = aqdp.u(idx,:);
b2 = aqdp.v(idx,:);
b3 = aqdp.w(idx,:);
t = aqdp.datenum(idx,:);

%correlations are low in upper bins. NaN these out.
ccrit = 70;
for i = badbins
    lowcor = aqdp.cor1(idx,i) <= ccrit;
    b1(lowcor,i) = NaN;
    lowcor = aqdp.cor2(idx,i) <= ccrit;
    b2(lowcor,i) = NaN;
    lowcor = aqdp.cor3(idx,i) <= ccrit;
    b3(lowcor,i) = NaN;
end

%rotate to the cross-shore direction
rot = heading*pi/180;
u = b1.*(ones(size(b1))*cos(rot)) + ...
         b2.*(ones(size(b2))*sin(rot));
v = -b1.*(ones(size(b1))*sin(rot)) + ...
        b2.*(ones(size(b2))*cos(rot));
w = b3;
rb = aqdp.rangebins;rbcm = h-rb*100;
sr = 8;
intv = 2; %minutes
avt = intv*sr*60;
ind = [1 avt:avt:length(w)];

for i = 1:length(ind)-1
    Uz = u(ind(i):ind(i+1),:);
    Uvar(i,:) = nanvar(Uz,0,1);
    Us(ind(i):ind(i+1),:) = Uz;
end


for i = 1:length(ind)-1
    Vz = v(ind(i):ind(i+1),:);
    Vvar(i,:) = nanvar(Vz,0,1);
    Vs(ind(i):ind(i+1),:) = Vz;
end

for i = 1:length(ind)-1
    Wz = w(ind(i):ind(i+1),:);
    Wvar(i,:) = nanvar(Wz,0,1);
    Ws(ind(i):ind(i+1),:) = Wz;
end

%convert to cm/s
Uvar = Uvar.*100;
Vvar = Vvar.*100;
Wvar = Wvar.*100;

c = linspecer(length(ind)-1,'sequential');
f1 = figure(1);
subplot(131)
for i = 1:length(ind)-1
    plot(Uvar(i,:),rbcm,'LineWidth',1.5,'Color',c(i,:),'marker','o')
    hold on
end
hold off
set(gca,'XLim',[0 10])
title('Along-Shore Variance')
xlabel('Velocity (cm/s)')
ylabel('Height Above Bed (cm)')

subplot(132)
for i = 1:length(ind)-1
    plot(Vvar(i,:),rbcm,'LineWidth',1.5,'Color',c(i,:),'marker','o')
    hold on
end
hold off
set(gca,'XLim',[0 10])
title('Cross-Shore Variance')
xlabel('Velocity (cm/s)')
ylabel('Height Above Bed (cm)')

subplot(133)
for i = 1:length(ind)-1
    plot(Wvar(i,:),rbcm,'LineWidth',1.5,'Color',c(i,:),'marker','o')
    hold on
end
hold off
set(gca,'XLim',[0 2])
title('Vertical Variance')
xlabel('Velocity (cm/s)')
ylabel('Height Above Bed (cm)')
suptitle('AD5116 - Above Bare-Bed')

figure(2)
subplot(131)
imagesc(t,rbcm,u')
set(gca,'Ydir','normal','YGrid','on')
datetickzoom('x','HH:MM','keepticks','keeplimits')
ylabel('Height Above Bed (cm)')
title('Along-Shore Velocities')

subplot(132)
imagesc(t,rbcm,v')
set(gca,'Ydir','normal','YGrid','on')
datetickzoom('x','HH:MM','keepticks','keeplimits')
ylabel('Height Above Bed (cm)')
title('Cross-Shore Velocities')

subplot(133)
imagesc(t,rbcm,w')
set(gca,'Ydir','normal','YGrid','on')
datetickzoom('x','HH:MM','keepticks','keeplimits')
ylabel('Height Above Bed (cm)')
title('Vertical Velocities')
suptitle('AD5116 - Above Bare-Bed')

clear Uvar Vvar Wvar aqdp
%%%%%%%

load('C:\Users\bkn5\Projects\Mekong_W2015\Data\Aquadopp\FSS\AD5117_9March2015.mat')
idx = find(aqdp.datenum >= start & aqdp.datenum <= stop);
h = 68.7; %inst hab (cm)
b1 = aqdp.u(idx,:);
b2 = aqdp.v(idx,:);
b3 = aqdp.w(idx,:);
t = aqdp.datenum(idx,:);

%correlations are low in upper bins. NaN these out.
ccrit = 70;
for i = 1:6
    lowcor = aqdp.cor1(idx,i) <= ccrit;
    b1(lowcor,i) = NaN;
    lowcor = aqdp.cor2(idx,i) <= ccrit;
    b2(lowcor,i) = NaN;
    lowcor = aqdp.cor3(idx,i) <= ccrit;
    b3(lowcor,i) = NaN;
end

%rotate to the cross-shore direction
rot = heading*pi/180;
u = b1.*(ones(size(b1))*cos(rot)) + ...
         b2.*(ones(size(b2))*sin(rot));
v = -b1.*(ones(size(b1))*sin(rot)) + ...
        b2.*(ones(size(b2))*cos(rot));
w = b3;
rb = aqdp.rangebins;rbcm = h-rb*100;
sr = 8;
avt = intv*sr*60;
ind = [1 avt:avt:length(w)];

for i = 1:length(ind)-1
    Uz = u(ind(i):ind(i+1),:);
    Uvar(i,:) = nanvar(Uz,0,1);
    Us(ind(i):ind(i+1),:) = Uz;
end

for i = 1:length(ind)-1
    Vz = v(ind(i):ind(i+1),:);
    Vvar(i,:) = nanvar(Vz,0,1);
    Vs(ind(i):ind(i+1),:) = Vz;
end

for i = 1:length(ind)-1
    Wz = w(ind(i):ind(i+1),:);
    Wvar(i,:) = nanvar(Wz,0,1);
    Ws(ind(i):ind(i+1),:) = Wz;
end
%convert to cm/s
Uvar = Uvar.*100;
Vvar = Vvar.*100;
Wvar = Wvar.*100;

phx = linspace(-10,10,10);
c = linspecer(length(ind)-1,'sequential');

f2 = figure(3);
subplot(131)
for i = 1:length(ind)-1
    plot(Uvar(i,:),rbcm,'LineWidth',1.5,'Color',c(i,:),'marker','o')
    hold on
end
plot(phx,ones(length(phx))*ph,'--k','LineWidth',1.5)
hold off
set(gca,'XLim',[0 10])
title('Along-Shore Variance')
xlabel('Velocity (cm/s)')
ylabel('Height Above Bed (cm)')

subplot(132)
for i = 1:length(ind)-1
    plot(Vvar(i,:),rbcm,'LineWidth',1.5,'Color',c(i,:),'marker','o')
    hold on
end
plot(phx,ones(length(phx))*ph,'--k','LineWidth',1.5)
hold off
set(gca,'XLim',[0 10])
title('Cross-Shore Variance')
xlabel('Velocity (cm/s)')
ylabel('Height Above Bed (cm)')

subplot(133)
for i = 1:length(ind)-1
    plot(Wvar(i,:),rbcm,'LineWidth',1.5,'Color',c(i,:),'marker','o')
    hold on
end
plot(phx,ones(length(phx))*ph,'--k','LineWidth',1.5)
hold off
set(gca,'XLim',[0 2])
title('Vertical Variance')
xlabel('Velocity (cm/s)')
ylabel('Height Above Bed (cm)')
suptitle('AD5117 - Above Pneumatophores')

figure(4)
subplot(131)
imagesc(t,rbcm,u')
set(gca,'Ydir','normal','YGrid','on')
datetickzoom('x','HH:MM','keepticks','keeplimits')
ylabel('Height Above Bed (cm)')
title('Along-Shore Velocities')

subplot(132)
imagesc(t,rbcm,v')
set(gca,'Ydir','normal','YGrid','on')
datetickzoom('x','HH:MM','keepticks','keeplimits')
ylabel('Height Above Bed (cm)')
title('Cross-Shore Velocities')

subplot(133)
imagesc(t,rbcm,w')
set(gca,'Ydir','normal','YGrid','on')
datetickzoom('x','HH:MM','keepticks','keeplimits')
ylabel('Height Above Bed (cm)')
title('Vertical Velocities')
suptitle('AD5117 - Above Pneumatophores')
