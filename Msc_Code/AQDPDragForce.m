%Calculate Drag Force from the seaward Aquadopp over the same period the
%VecPros were recording

%need to plot Northing/Easting component vectors of velocity from the
%mudflat Aquadopp to determine the direct onshore velocity (u) to calculate
%the relationship between average drag (Fd) and TKE dissipation (epsilon),
%given in Dalrymple et al., 1989

%This definition is:
%epsilon_d = F_d*u = (1/2)*rho*a*C_d*|u^3| where F_d*u and |u^3| are time
%averaged WPF*velocity and velocity cubed, respectively. Definitions of C_d
%and a are C_d = 1 and a = 0.5m^-1.
% clear
load HR3_12March2015_f_pad.mat
dname = 'FSS33'; %deployment
name = 'ADHR3'; %instrument
% savefigdir = 'c:\Users\bkn5\Projects\Mekong_W2015\Figures\Aquadopps\HR3_9March2015\';
intv = 10; %minutes
start = datenum(2015,03,10,14,45,00);
stop = datenum(2015,03,10,16,25,00);
%NOTE: Aquadopps recorded in local (ICT) timezone
% plot(aqdp.datenum,aqdp.pressure)
% datetickzoom('x','dd HH:MM:SS','keepticks','keeplimits')
savedatfile = 1;
%Initialize Global Variables %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ind = find(aqdp.datenum >= start & aqdp.datenum <=stop);
u = aqdp.u(ind,:);
v = aqdp.v(ind,:);
Cd = 1; %drag coef
a = 0.5; %m^-1

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
p = aqdp.pressure(ind);
t = aqdp.temperature(ind);
s = repmat(20,length(t),1); %assume constant salinity of 20PSU
rho = SeaDensity(s,t,p);
rho = nanmean(rho);
%Fill gaps in timeseries if they exist %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if cmgidgaps(u) > 0
    disp(['Found ' num2str(cmgidgaps(u)) ' gaps in u time-series'])
    nlin = 1E2;
    maxgaps = 1E5;
    u = cmgbridge(u,nlin,maxgaps,maxgaps);
end
if cmgidgaps(v) > 0
    disp(['Found ' num2str(cmgidgaps(v)) ' gaps in v time-series'])
    v = cmgbridge(v,nlin,maxgaps,maxgaps);
end
if cmgidgaps(u) > 0 || cmgidgaps(v) > 0
    u(isnan(u)) = 0;
    v(isnan(v)) = 0;
    disp(['Number of gaps in velocity t-s remaining: ' num2str(cmgidgaps(u))])
else
    disp('Gaps filled')
end

%Compute depth averages %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[m,n] = size(u);
ud = zeros(m,1);vd = zeros(m,1);
for i = 1:m
    ud(i,1) = nanmean(u(i,1:n));
    vd(i,1) = nanmean(v(i,1:n));
end

%Calculate regression %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
mm = min(ud);
ma = max(ud);
xes = linspace(mm,ma,length(ud));
pf = polyfit(ud,vd,1);
pv = polyval(pf,xes);
theta = tan(pf(1)); %angle of slope is the tangent
deg = rad2deg(theta);
%0 deg is due W/E; subtract 90deg from this value to get due north, then
%adjust for the x-shore transect angle
if deg < 0
    xsdeg = 270+abs(deg);
    xsdeg = 350-xsdeg;
elseif deg > 0
    xsdeg = 270-deg;
    xsdeg = 350-xsdeg;
end

disp(['Value to rotate u,v velocity data by: ' num2str(xsdeg)])

%lets try it
[u,v] = mag2truenorth(u,v,xsdeg);
%The mag2truenorth script rotates cw, and the resultant v component
%(north/south) velocities are onshore/offshore.

%Re-compute depth averages %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ud2 = zeros(m,1);vd2 = zeros(m,1);
for i = 1:m
    ud2(i,1) = nanmean(u(i,1:n));
    vd2(i,1) = nanmean(v(i,1:n));
end

%Rotate line %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
mm = min(ud2);
ma = max(ud2);
xes = linspace(mm,ma,length(ud));
xsrad = xsdeg*pi/180;
[th,r] = cart2pol(xes,pv);
th = th-xsrad;
[xes2,pv2] = pol2cart(th,r);
%double check angle
pf = polyfit(xes2,pv2,1);
theta = tan(1/(pf(1))); %angle of slope is the tangent
deg = rad2deg(theta);
if deg < 0
    rotdeg = 360+deg;
elseif deg > 0 
    rotdeg = 360-deg;
end
disp(['Onshore flow direction is ' num2str(rotdeg) ' degrees'])

%Plot results %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
f1 = figure;
set(f1,'PaperOrientation','portrait',...
    'position',[400 200   600   400]);
set(gcf,'color','w','PaperPositionMode','auto')
h(1) = plot(ud,vd,'or');
hold on
h(2) = plot(xes,pv,'LineWidth',1.5,'Color','k');
h(3) = plot(ud2,vd2,'ob');
h(4) = plot(xes2,pv2,'LineWidth',1.5,'Color','k');
text(0.05,0.05,['Rotated angle: ' sprintf('%c', char(176))...
    sprintf('%0.2f',rotdeg) sprintf('%c', char(176))],...
    'units','normalized')
xlabel('\bf\itEast Velocity (m/s)')
ylabel('\bf\itNorth Velocity (m/s)')
title(['\bf\itDepth Averaged Velocities, ' name ' ' dname])
axis equal
legend([h(1),h(3)],{'Unrotated';'Rotated'},'location','northeast')
% prompt = 'Save Figure? [y/n] ';
% result = input(prompt,'s');
% if strcmp(result,'y') || strcmp(result,'yes');
%     fpath = savefigdir;fname = [name '_uvdavel'];
%     export_fig([fpath fname],'-png','-m1','-r900','-opengl')
%     disp(['Figure ' fname '.png saved'])
% end
% if strcmp(result,'n') || strcmp(result,'no');
% end

DragF = struct();
fs = str2double(aqdp.metadata.samprate(1:2)); %sampling frequency
avt = fs*intv*60; %# samples in interval
idx = avt:avt:length(vd2);
idxx = [1 idx];
Fd = zeros(1,length(idx));
u = zeros(length(idx),1);
time = zeros(1,length(idx));
it = 1;
while it < length(idxx)
    ind = idxx(it):idxx(it+1);
    %Calculate equivalent drag force %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Definition: Epsilon_d = F_d*v = (1/2)*rho*a*Cd*|v^3|
    v = rms(vd2(ind,:));
    Fd(:,it) = (1/2)*rho*a*Cd.*v.^3;
    u(it,:) = v;
    time(:,it) = intv*it; 
    it = it+1; 
end
DragF.DepStart = [datestr(start,'dd-mm-yyyy HH:MM:SS') ' ICT'];
DragF.DepStop = [datestr(stop,'dd-mm-yyyy HH:MM:SS') ' ICT'];
DragF.avgint = [num2str(intv) ' minute averages'];
DragF.cmt = 'Each column represents a new interval';
DragF.u = u;
DragF.rho = rho;
DragF.Cd = Cd;
DragF.a = a;
DragF.Fd = Fd;
DragF.Time = time;
if savedatfile
    %save file
    filename = [name '_' dname];
    sfname = ['Fd' filename];
    save(sfname,'DragF','-v7.3')
    disp(['Drag force data file saved as ' sfname])
end