close all
clear

pw = 'D:\Projects\Mekong_W2015\Data\Aquadopp\DPS2\';
files = dir([pw '*15March2015.mat']);files = {files.name};
names = {'WU'};
symb = {'d'};
scale = 0.01;

load([pw files{4}])
disp(['Loading ' files{4}])
lon = aqdp.metadata.lon;
lat = aqdp.metadata.lat;
u = nanmean(aqdp.u,2);u(isnan(u)) = 0;
v = nanmean(aqdp.v,2);v(isnan(v)) = 0;
time = aqdp.datenum;
pres = cmgbridge(aqdp.pressure,100,1000,10000);

figure(1)
n = length(u);
pd(1) = plot(lon,lat,'Marker',symb{1},...
    'MarkerSize',8,'LineWidth',1,'Color','k'); hold on

c = jet(n);
% posX = cumsum([lon ; u(:).*scale]);
% posY = cumsum([lat ; v(:).*scale]);
for i = 1:n-1
    quiver(lon,lat,v(i),u(i),scale,...
        'color',c(i,:))
    
    
    
%     line(posX(i:i+1),posY(i:i+1),...
%         'LineWidth',1.5,...
%         'Color',c(i,:));
end

savefigdir = 'd:\Projects\Mekong_W2015\Figures\Environment\';
xlabel('lon')
ylabel('lat')
title('Progressive Vectors')
legend(pd,names)
export_fig([savefigdir 'Aarons_PVD'],'-pdf','-nocrop')
