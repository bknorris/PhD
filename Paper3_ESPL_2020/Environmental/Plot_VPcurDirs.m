figure(1)


u = rxy(25000:25250,1);
v = rxy(25000:25250,2);

water = [186 213 247]/255;
mud = [214 199 174]/255;
veg = [195 204 169]/255;
%Plot fringe line on figures
fdir = 'e:\Mekong_W2015\ArcGis\GeoRef\';
fdat = fopen([fdir 'SW_fringe_vf.csv']);
fgetl(fdat);
fdat = textscan(fdat,'%n%n%n','delimiter',',');
flon = fdat{:,1};flat = fdat{:,2};
subplot(131)
area(flon,flat,9.45,...
    'LineWidth',1.5,...
    'facecolor',water,...
    'edgecolor',mud), hold on
area(flon,flat,9.6,...
    'LineWidth',1.5,...
    'facecolor',veg,...
    'edgecolor',mud)
plot(flon,flat,':',...
    'LineWidth',1.5,...
    'Color','k')
symb = {'d'};
scale = 0.01;
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
subplot(132)
for i = 1:n-1
    plot(vpt(i),u(i),'.','color',c(i,:)),hold on
end
subplot(133)
for i = 1:n-1
    plot(vpt(i),v(i),'.','color',c(i,:)),hold on
end
