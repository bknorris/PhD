%Figure out the current rotations for the DPS (2015) experiment

close all
clear

pw = 'D:\Projects\Mekong_W2015\Data\Aquadopp\DPS2\';
files = dir([pw '*15March2015.mat']);files = {files.name};
names = {'5116';'5117';'HR3';'UW';'V5109';'VP3'};
symb = {'o';'^';'s';'d'};
scale = 0.01;
start = datenum(2015,03,13,16,20,00);
stop = datenum(2015,03,13,21,43,00);
for ff = 1:length(files)
    load([pw files{ff}])
    disp(['Loading ' files{ff}])
    
    ind = find(aqdp.datenum >= start & aqdp.datenum <= stop);
    if ff < 4
        u = cmgbridge(nanmean(aqdp.u(ind,:),2),100,1000,10000);
        v = cmgbridge(nanmean(aqdp.v(ind,:),2),100,1000,10000);
    else
        u = cmgbridge(nanmean(aqdp.u(ind,:),2),100,1000,10000);
        v = cmgbridge(nanmean(aqdp.v(ind,:),2),100,1000,10000);
    end
    lat = aqdp.metadata.lat;
    lon = aqdp.metadata.lon;
    time = aqdp.datenum(ind);
    pres = cmgbridge(aqdp.pressure(ind),100,1000,10000);
    fs = 8;
    intv = 7.5;
    avt  = fs*intv*60;
    id = [1 avt:avt:length(u)];
    
    figure(1)
    if ff < 4
        c = jet(length(id));
        for i = 1:length(id)-1
            quiver(lon,lat,mean(u(id(i):id(i+1))),mean(v(id(i):id(i+1))),scale,...
                'Color',c(i,:));hold on
        end
    else
        c = jet(length(u));
        for i = 1:length(u)
            quiver(lon,lat,u(i),v(i),scale,...
                'Color',c(i,:));hold on
        end
    end
    text(lon,lat,names{ff},'FontSize',10)
    
    figure(2)
    if ff < 4
        c = jet(length(id));
        for i = 1:length(id)-1
            pp(ff) = plot(time(id(i)),mean(pres(id(i):id(i+1))),...
                'Color','k','MarkerFaceColor',c(i,:),...
                'Marker',symb{ff});hold on
        end
    else
        c = jet(length(u));
        for i = 1:length(u)
            pp(ff) = plot(time(i),pres(i),...
                'Color','k','MarkerFaceColor',c(i,:),...
                'Marker',symb{ff});hold on
        end
    end
    datetickzoom('x','HH:MM:SS','keepticks','keeplimits')
    
    figure(3)
    n = length(u);
    pd(ff) = plot(lon,lat,'Marker',symb{ff},...
        'MarkerSize',8,'LineWidth',1,'Color','k'); hold on
    if ff < 4
        c = jet(length(id));
        for i = 1:length(id)-1
            umean(i) = mean(u(id(i):id(i+1)));
            vmean(i) = mean(v(id(i):id(i+1)));
            posX = cumsum([lon ; umean(:).*scale]);
            posY = cumsum([lat ; vmean(:).*scale]);
        end
        for i = 1:length(id)-1
            line(posX(i:i+1),posY(i:i+1),...
                'LineWidth',1.5,...
                'Color',c(i,:));
        end
    else
        u(isnan(u)) = 0;v(isnan(v)) = 0;
        c = jet(n);
        posX = cumsum([lon ; u(:).*scale]);
        posY = cumsum([lat ; v(:).*scale]);
        for i = 1:n-1
            line(posX(i:i+1),posY(i:i+1),...
                'LineWidth',1.5,...
                'Color',c(i,:));
        end
    end
end

%Load the Vector
pw = 'D:\Projects\Mekong_W2015\Data\Vector\DPS2\';
files = dir([pw '*130315.mat']);files = {files.name};
symb = {'p'};

load([pw files{1}])
disp(['Loading ' files{1}])

ind = find(ADV.datetime >= start & ADV.datetime <= stop);
u = cmgbridge(nanmean(ADV.U(ind),2),100,1000,10000);
v = cmgbridge(nanmean(ADV.V(ind),2),100,1000,10000);
lat = ADV.Metadata.inst_lat;
lon = ADV.Metadata.inst_lon;
time = ADV.datetime(ind);
pres = cmgbridge(ADV.Pres(ind),100,1000,10000);
fs = 32;
intv = 8;
avt  = fs*intv*60;
id = [1 avt:avt:length(u)];

f1 = figure(1);
c = jet(length(id));
for i = 1:length(id)-1
    quiver(lon,lat,mean(u(id(i):id(i+1))),mean(v(id(i):id(i+1))),scale,...
        'Color',c(i,:));hold on
end

text(lon,lat,names{5},'FontSize',10)

f2 = figure(2);
c = jet(length(id));
for i = 1:length(id)-1
    pp(5) = plot(time(id(i)),mean(pres(id(i):id(i+1))),...
        'Color','k','MarkerFaceColor',c(i,:),...
        'Marker',symb{1});hold on
end

%     datetickzoom('x','HH:MM:SS','keepticks','keeplimits')

f3 = figure(3);
n = length(u);
pd(5) = plot(lon,lat,'Marker',symb{1},...
    'MarkerSize',8,'LineWidth',1,'Color','k'); hold on
c = jet(length(id));
for i = 1:length(id)-1
    umean(i) = mean(u(id(i):id(i+1)));
    vmean(i) = mean(v(id(i):id(i+1)));
    posX = cumsum([lon ; umean(:).*scale]);
    posY = cumsum([lat ; vmean(:).*scale]);
end
for i = 1:length(id)-1
    line(posX(i:i+1),posY(i:i+1),...
        'LineWidth',1.5,...
        'Color',c(i,:));
end

%Load the Vectrino (VP3) - above the canopy
pw = 'D:\Projects\Mekong_W2015\Data\Vectrino\13March2015\';
file = 'VP3_xyzvels.mat';
load([pw file])
symb = {'h'};
ind = find(VP3.datetime >= start & VP3.datetime <= stop);

%rotate the VP to north X is north-south, Y is east-west
heading = 22; %rotate CCW
rot = heading*pi/180;
v = VP3.x.*(ones(size(VP3.x))*cos(rot)) + ...
    VP3.y.*(ones(size(VP3.y))*sin(rot));
u = -VP3.y.*(ones(size(VP3.y))*sin(rot)) + ...
    VP3.x.*(ones(size(VP3.x))*cos(rot));


u = cmgbridge(nanmean(u(ind),2),100,1000,10000);
v = cmgbridge(nanmean(v(ind),2),100,1000,10000);
lat = 9.565494;
lon = 106.292603;
fs = 50;
intv = 8;
avt  = fs*intv*60;
id = [1 avt:avt:length(u)];

f1 = figure(1);
c = jet(length(id));
for i = 1:length(id)-1
    quiver(lon,lat,mean(u(id(i):id(i+1))),mean(v(id(i):id(i+1))),scale,...
        'Color',c(i,:));hold on
end
text(lon,lat,names{6},'FontSize',10)

f3 = figure(3);
n = length(u);
pd(6) = plot(lon,lat,'Marker',symb{1},...
    'MarkerSize',8,'LineWidth',1,'Color','k'); hold on
c = jet(length(id));
for i = 1:length(id)-1
    umean(i) = mean(u(id(i):id(i+1)));
    vmean(i) = mean(v(id(i):id(i+1)));
    posX = cumsum([lon ; umean(:).*scale]);
    posY = cumsum([lat ; vmean(:).*scale]);
end
for i = 1:length(id)-1
    line(posX(i:i+1),posY(i:i+1),...
        'LineWidth',1.5,...
        'Color',c(i,:));
end

savefigdir = 'd:\Projects\Mekong_W2015\Figures\Environment\';
figure(1)
xlabel('lon')
ylabel('lat')
title('Current Vectors, DPS-2')
% export_fig([savefigdir 'NW_allCurvec'],'-pdf','-nocrop')
figure(2)
xlabel('Time on 13/03/2015')
ylabel('dbar')
title('Pressure Records, DPS-2 Multiple Instruments')
legend(pp,names)
% export_fig([savefigdir 'NE_allPres'],'-pdf','-nocrop')
figure(3)
xlabel('lon')
ylabel('lat')
title('Progressive Vectors, DPS-2 Multiple Instruments')
legend(pd,names)
% export_fig([savefigdir 'NE_allPVDs'],'-pdf','-nocrop')




