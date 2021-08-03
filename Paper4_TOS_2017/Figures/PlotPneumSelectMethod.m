clear
folder = 'd:\Projects\Mekong_W2015\DataAnalysis\Paper2\VegStats';
% folder = uigetdir;
fname = dir([folder '\' '*stats_fixed.mat']);fname = {fname.name};
vegdat = struct();
savefigdir = 'd:\Projects\Mekong_W2015\Figures\Paper2\Phi-Eps\';

%Load the VP positions data
vpdir = 'd:\Projects\Documents\Writing\DataReports\';
load([vpdir 'VP_positions.mat'])
VPpos.Qname = VPpos.Qname(22:end);
VPpos.VP_X = VPpos.VP_X(22:end);
VPpos.VP_Y = VPpos.VP_Y(22:end);
rotx = [0 180 180 0 0 180 0 180];
roty = [0 0 0 0 180 0 0 180];
rotz = [90 0 0 270 90 0 -90 0];
k = 1;

%%%Load the Data
disp(['Loading ' folder '\' fname{k}])
load([folder '\' fname{k}])
qname = regexprep(fname{k},'stats_fixed.mat','');
DATA = rmfield(DATA,'info'); %remove info field for processing

%%%Load the correct VP data and scatterplot
vpid = zeros(length(VPpos.Qname),1);
for b = 1:length(VPpos.Qname)
    vpid(b) = strcmp(VPpos.Qname{b},qname);
end
vpx = (VPpos.VP_X.*vpid)';vpx(vpx == 0) = [];
vpy = (VPpos.VP_Y.*vpid)';vpy(vpy == 0) = [];

%%%First, plot the entire quadrat as an overview
f1 = figure(1);
set(f1,'PaperOrientation','portrait',...
    'position',[200 300   800   600]);
set(gcf,'color','w','PaperPositionMode','auto')
%Create rotation matrix, rotate about z-axis if specified, flip
%axes if specified
thx = rotx(k)*pi/180;thy = roty(k)*pi/180;thz = rotz(k)*pi/180;
Tx = [cos(2*thx) sin(2*thx); sin(2*thx) -cos(2*thx)];
Ty = [-cos(2*thy) sin(2*thy); sin(2*thy) cos(2*thy)];
Tz = [cos(thz) -sin(thz); sin(thz) cos(thz)];

%Adjust VP positions based on rotations and reflections
o = repmat([0.5;0.5],1,length(vpx)); %origin
vpxy = [vpx;vpy];
if rotz(k) ~= 0
    vpxy = Tz*(vpxy-o)+o;
end
if rotx(k) ~= 0
    vpxy = Tx*(vpxy-o)+o;
end
if roty(k) ~= 0
    vpxy = Ty*(vpxy-o)+o;
end

l = length(DATA.layer1.Rcenter);
x = DATA.layer1.Rcenter(:,1)';
y = DATA.layer1.Rcenter(:,2)';

%Rotate & Reflect coordinates
o = repmat([0.5;0.5],1,length(x));
xy = [x;y];
if rotz(k) ~= 0
    xy = Tz*(xy-o)+o;
end
if rotx(k) ~= 0
    xy = Tx*(xy-o)+o;
end
if roty(k) ~= 0
    xy = Ty*(xy-o)+o;
end

for i = 1:l
    xc = xy(1,i);
    yc = xy(2,i);
    r = mean(DATA.layer1.Rradii(i,:));
    ang=0:0.01:2*pi;
    xp=r*cos(ang);
    yp=r*sin(ang);
    p = plot(xc+xp,yc+yp);set(p,'linewidth',2,...
        'Color','k')
    hold on
end
axis equal
grid on
set(gca,'Xlim',[0 1],'Ylim',[0 1],'grid',':',...
    'XTick',0:0.1:1,'YTick',0:0.1:1,...
'LineWidth',1.5,'FontSize',13,'FontName','Arial',...
'Color',[0.9 0.9 0.9])
xlabel('Easting (m)','FontSize',14),ylabel('Northing (m)','FontSize',14)
title(['Fall 2014 - ' qname],'FontSize',14)

%%%Draw a rectangle around the desired points (for help: 'help rbbox')
hold on
c = brewermap(4,'YlGnBu');
colormap(c)
rectangle('position',[0.5,0.44,0.2,0.2],'EdgeColor',c(1,:),'LineWidth',2)
rectangle('position',[0.52,0.43,0.2,0.2],'EdgeColor',c(2,:),'LineWidth',2)
rectangle('position',[0.57,0.44,0.2,0.2],'EdgeColor',c(3,:),'LineWidth',2)
rectangle('position',[0.54,0.45,0.2,0.2],'EdgeColor',c(4,:),'LineWidth',2)
xlo = 0.54;xhi = 0.74;
ylo = 0.45;yhi = 0.65;
%Plot VP Positions
plot(vpxy(1,:),vpxy(2,:),'o',...
    'Color','r',...
    'MarkerSize',18,...
    'LineWidth',1.5)
plot(vpxy(1,:),vpxy(2,:),'^',...
    'MarkerSize',15,...
    'MarkerEdgeColor','k',...
    'MarkerFaceColor','k',...
    'LineWidth',1.5)

% export_fig([savefigdir 'PneumSelectLg'],'-pdf','-nocrop')

%%%Crop DATA layer1 based on answer inputs, then replot
xid = find(xy(1,:) >= xlo & xy(1,:) <= xhi);
yid = find(xy(2,:) >= ylo & xy(2,:) <= yhi);
id = intersect(xid,yid);


f2 = figure(2);
set(f2,'PaperOrientation','portrait',...
    'position',[200 300   800   600]);
set(gcf,'color','w','PaperPositionMode','auto')

%%%Plot VP Positions
%Plot VP Positions
plot(vpxy(1,:),vpxy(2,:),'o',...
    'Color','r',...
    'MarkerSize',18,...
    'LineWidth',1.5), hold on
plot(vpxy(1,:),vpxy(2,:),'^',...
    'MarkerSize',15,...
    'MarkerEdgeColor','k',...
    'MarkerFaceColor','k',...
    'LineWidth',1.5)

l = length(id);
R = zeros(l,1);
X = zeros(l,1);
Y = zeros(l,1);
for i = 1:l
    xc = xy(1,id(i));
    yc = xy(2,id(i));
    r = mean(DATA.layer1.Rradii(id(i),:));
    ang=0:0.01:2*pi;
    xp=r*cos(ang);
    yp=r*sin(ang);
    p = plot(xc+xp,yc+yp);set(p,'linewidth',2,...
        'Color','k')
    hold on
    R(i,:) = r;
    X(i,:) = xc;
    Y(i,:) = yc;
end
axis equal
grid on
set(gca,'Xlim',[xlo xhi],'Ylim',[ylo yhi],...
    'XTick',xlo:0.1:xhi,...
    'YTick',ylo:0.1:yhi,...
    'grid',':',...
    'LineWidth',1.5,...
    'Color',[0.9 0.9 0.9],...
    'FontSize',13,'FontName','Arial')
xlabel('Easting (m)','FontSize',14),ylabel('Northing (m)','FontSize',14)
title(['Fall 2014 - ' qname ' Cropped'],'FontSize',14)
hold off

% export_fig([savefigdir 'PneumSelectSm'],'-pdf','-nocrop')

