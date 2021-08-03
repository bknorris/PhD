%Make plots of current vectors to determine the dominant angle of onshore
%flow throughout the experiment. First calculate depth averaged velocities,
%then plot mean velocity pairs (X,Y) corresponding to a given time window.
%Will need to convert beam coordinates to XYZ.
close all
clear
load VP3_090315.mat 
dname = 'FSS3'; %deployment
name = 'VP3'; %instrument
iname = VPRO.Config.Ins_name;
savefigdir = 'c:\Users\bkn5\Projects\Mekong_W2015\Figures\Vectrinos\Mar9\VP3\'; 

%convert to XYZ coordinates
[VPRO.Data,VPRO.Config] = VPro_coordinateTransform(VPRO.Data,VPRO.Config,'bx');
sr = VPRO.Config.sampleRate;

%compute depth averages
[m,n] = size(VPRO.Data.Profiles_VelX);
x = zeros(m,1);y = zeros(m,1);z1 = zeros(m,1);z2 = zeros(m,1);
for i = 1:m
    x(i,1) = nanmean(VPRO.Data.Profiles_VelX(i,1:n));
    y(i,1) = nanmean(VPRO.Data.Profiles_VelY(i,1:n));
%     z1(i,1) = nanmean(VPRO.Data.Profiles_VelZ1(i,1:n));
%     z2(i,1) = nanmean(VPRO.Data.Profiles_VelZ2(i,1:n));
end
clearvars VPRO

%window
intv = 0.5; %averaging window in minutes
avt = 60*intv*sr;
ind = avt:avt:m;    
ind = [1 ind m]; %include remainder
p = length(ind)-1;
xav = zeros(p,1);yav = zeros(p,1);%z1av = zeros(p,1);z2av = zeros(p,1);
for i = 1:p
    idx = ind(i):ind(i+1);
    xav(i,1) = nanmean(x(idx));
    yav(i,1) = nanmean(y(idx));
%     z1av(i,1) = nanmean(z1(idx));
%     z2av(i,1) = nanmean(z2(idx));
end
% zav = (z1av+z2av)./2;

%calculate regression
pf = polyfit(xav,yav,1);
pv = polyval(pf,xav);
theta = tan(pf(1)); %angle of slope is the tangent
degoffn = rad2deg(theta); %degrees from north
if degoffn > 0 
    degoffn = 360 - degoffn;
elseif degoffn < 0
    degoffn = abs(degoffn);
end

%plot
f1 = figure;
set(f1,'PaperOrientation','portrait',...
    'position',[400 200   600   400]);
set(gcf,'color','w','PaperPositionMode','auto')
h(1) = plot(xav,yav,'or');
hold on
h(2) = plot(xav,pv,'LineWidth',1.5,'Color','k');
text(0.1,0.1,['Degrees off Mag N ' sprintf('%0.2f',degoffn) sprintf('%c', char(176))],...
    'units','normalized')
xlabel(['\bf\it' '$\bar{x}$' ' (m/s)'],'interpreter','latex')
ylabel(['\bf\it' '$\bar{y}$' ' (m/s)'],'interpreter','latex')
title(['\bf\itDepth Averaged Velocities, ' iname ' ' dname ' in ' ...
    num2str(intv) ' minute averages'])
axis equal
prompt = 'Save Figure? [y/n] ';
result = input(prompt,'s');
if strcmp(result,'y') || strcmp(result,'yes');
    fpath = savefigdir;fname = [name '_XYdavel'];
    export_fig([fpath fname],'-png','-m1','-r900','-opengl')
    disp(['Figure ' fname '.png saved'])
end
if strcmp(result,'n') || strcmp(result,'no');
end

