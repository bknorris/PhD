%using the difference map created in ArcGIS, load the raster .ascii file,
%constrain the x and y coordinates to a unit vector limits [0 1] meters and
%plot.
clear
dirc = 'D:\Projects\Mekong_W2015\Images\Quadrats\GIS';
qdat = 'D:\Projects\Mekong_W2015\Images\Quadrats\Q2\Reconstructions\2A\Analysis\';
fdir = 'D:\Projects\Mekong_W2015\Images\Quadrats\Q2\';
fname = 'twoatwobdiff.txt';cd(dirc)
name = 'Q2A2BdiffMap';
[z,r] = arcgridread(fname);
%measure one dimension of the cloud on cloud compare. This value is equal
%to 1m, or, to convert elevations to meters, divide Z by this value.
ccscale = 6.657;
Z = (z./ccscale).*100; %*100 converts to cm
Z = rot90(Z,-1); %rotate Z so it is in the same view as the orthophotos

%create unit vector axis
[m,n] = size(Z);
x = linspace(0,1,n);X = repmat(x,m,1);
y = linspace(1,0,m);Y = repmat(y,n,1)';


%load quadrat data
load([qdat 'Q2Astats_fixed.mat'])

%contour using contourf
f1 = figure(1);
set(f1,'PaperOrientation','portrait',...
    'position',[400 200   800   800]);
set(gcf,'color','w','PaperPositionMode','auto')

[c,h] = contourf(X,Y,Z);
set(h,'EdgeColor','none')

hold on
%calculate radii of pneumatophores at level 1 of the reconstruction
nn = length(DATA.layer1.Rcenter);
for j = 1:nn %loop through rows for the Radius and Center points of each layer
    cx = DATA.layer1.Rcenter(j,1);
    cy = DATA.layer1.Rcenter(j,2);
    
    %transformation matrix, rotate 90deg CCW
    R = [0 -1;1 0];
    rot = R*[cx;cy];
    cx = abs(rot(1));cy = rot(2);
    cr = mean(DATA.layer1.Rradii(j,:));
    ang=0:0.01:2*pi;
    xp=cr.*cos(ang);
    yp=cr.*sin(ang);
    %plot circles with the same radius as the vegetation
    H = patch(cx+xp,cy+yp,1);
    set(H,'FaceColor',[1 0 0],'EdgeColor','none')
end
hold off
set(gca,'LineWidth',1.5,'FontSize',20)
xlabel('\bf\itAlong Shore Distance (m)')
ylabel('\bf\itCross-Shore Distance (m)')
caxis([-5 5])
cb = colorbar;set(cb,'LineWidth',1.5,'FontSize',20)
ylabel(cb,'\bf\itBed Level Change (cm)')
% export_fig([fdir name],'-pdf','-m1','-r900')