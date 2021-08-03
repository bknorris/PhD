%Transform Vector coordinates in BEAM to ENU using some assumed values for
%heading, pitch and roll.

%This scipt will do the rotations, then plot the results. This is only for
%checking to see if the data are oriented in space in the way that you
%expect, and should not be used for bulk data processing or analysis (due
%to the assumptions made).

clear
close all
%load a concatenated (or raw) VecPro file
fname = 'VP3_05032015FC.mat';
israwfile = 0; %set this flag to 1 if the file being loaded is a raw .mat file
hcaxis = 0.25; %horizontal color axis for plots (U and V)
vcaxis = 0.025; %vertical color axis for plots (W)
vname = 'VP3'; %vectrino name, either VP1, VP2, or VP3
savefigdir = 'c:\Users\bkn5\Projects\Mekong_W2015\Figures\Vectrinos\Mar5\VP3\FlumeControl\10minavs\'; 
%fig directory
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

load(fname)
if israwfile
    VPRO.Data = Data;
    VPRO.Data.Time = Data.Profiles_HostTimeMatlab;
    VPRO.Config = Config;
    clearvars Data Config
end
%all VecPro files we recorded have velocities in BEAM coordinates
disp('Transforming VecPro velocities in BEAM coordinates to ENU')

if isfield(VPRO.Config,'nCells')
    nCells = VPRO.Config.nCells;
    nCalibratedCells = size(VPRO.Config.ProbeCalibration_calibrationMatrix,1);
    originalPrefix = '';
elseif isfield(VPRO.Config,'Original_nCells')
    nCells = VPRO.Config.Original_nCells;
    nCalibratedCells = size(VPRO.Config.ProbeCalibration_calibrationMatrix,1);
    originalPrefix = 'Original_';
end

%preallocate variables
VPRO.Data.Profiles_U = NaN*zeros(size(VPRO.Data.Profiles_VelBeam1));
VPRO.Data.Profiles_V = NaN*zeros(size(VPRO.Data.Profiles_VelBeam1));
VPRO.Data.Profiles_W = NaN*zeros(size(VPRO.Data.Profiles_VelBeam1));

%average the Z beams together
Z = (VPRO.Data.Profiles_VelBeam3+VPRO.Data.Profiles_VelBeam4)/2;

%define the rotation parameters
heading = 0; %magnetic north
pitch = 0;
roll = 0;

%now do the rotations
for cell = 1:nCells
    if cell <= nCalibratedCells
        T = reshape(VPRO.Config.([originalPrefix 'ProbeCalibration_calibrationMatrix'])(cell,:),4,4)';
    else
        T = eye(4);
    end
    V = double([VPRO.Data.Profiles_VelBeam1(:,cell)';
        VPRO.Data.Profiles_VelBeam2(:,cell)';
        Z(:,cell)']);
    T = T(1:3,1:3); %reduce the transformation matrix to 3x4
    hh = pi*(heading-90)/180;
    pp = pi*pitch/180;
    rr = pi*roll/180;
    H = [cos(hh) sin(hh) 0; -sin(hh) cos(hh) 0; 0 0 1];
    
    % Make tilt matrix
    P = [cos(pp) -sin(pp)*sin(rr) -cos(rr)*sin(pp);...
        0             cos(rr)          -sin(rr);  ...
        sin(pp) sin(rr)*cos(pp)  cos(pp)*cos(rr)];
    % Make resulting transformation matrix depending on coord system
    R = H*P*T;
    V = R*V;
    VPRO.Data.Profiles_U(:,cell) = V(1,:)';
    VPRO.Data.Profiles_V(:,cell) = V(2,:)';
    VPRO.Data.Profiles_W(:,cell) = V(3,:)';
end
VPRO.Config.([originalPrefix 'coordSystem']) = 'ENU';
VPRO.Data = rmfield(VPRO.Data,{'Profiles_VelBeam1','Profiles_VelBeam2','Profiles_VelBeam3','Profiles_VelBeam4'});
disp('Transformed BEAM to ENU')
%plot dat shit
time = VPRO.Data.Profiles_HostTimeMatlab;
b1 = VPRO.Data.Profiles_U;
b2 = VPRO.Data.Profiles_V;
b3 = VPRO.Data.Profiles_W;
titles = {'U';'V';'W'};

rangebins = VPRO.Data.Profiles_Range;
date = datestr(time(1),'dd/mm/yyyy');
max1 = hcaxis;
min1 = -hcaxis;
step = (-min1+max1)/4;
max2 = vcaxis;
min2 = -vcaxis;
step2 = (-min2+max2)/4;

f1 = figure;
set(f1,'PaperOrientation','portrait',...
    'position',[200 200   1250   550]);
ax = tight_subplot(3,1,0.05,[0.1 0.05],[0.1 0.05]);
axes(ax(1))
imagesc(time,flipud(rangebins),flipud(b1'))
title(['\bf\it' titles{1} ' Velocities'])
caxis([min1 max1])
c = colorbar;
cbh = get(c,'Title');
titleString = '[m/s]';
set(cbh ,'String',titleString);
set(c,'YTick',[min1:step:max1]);
g = gca;
set(g,'XTickLabel',[])
ylabel('\bf\itrange (m)')

axes(ax(2))
imagesc(time,flipud(rangebins),flipud(b2'))
title(['\bf\it' titles{2} ' Velocities'])
caxis([min1 max1])
c = colorbar;
cbh = get(c,'Title');
titleString = '[m/s]';
set(cbh ,'String',titleString);
set(c,'YTick',[min1:step:max1]);
g = gca;
set(g,'XTickLabel',[])
ylabel('\bf\itrange (m)')

axes(ax(3))
imagesc(time,flipud(rangebins),flipud(b3'))
title(['\bf\it' titles{3} ' Velocities'])
caxis([min2 max2])
c = colorbar;
cbh = get(c,'Title');
titleString = '[m/s]';
set(cbh ,'String',titleString);
set(c,'YTick',[min2:step2:max2]); %#ok<*NBRAK>
g = gca;
set(g,'XTickLabel',[])
ylabel('\bf\itrange (m)')
xlabel(['\bf\itTime on ' date])
set(ax,...
    'YTick',[0.04:0.01:0.07],...
    'YLim', [0.04 0.07],...
    'XGrid','on',...
    'GridLineStyle','-',...
    'YDir','normal',...
    'TickDir','in',...
    'TickLength',[.01 .01],...
    'XMinorTick','off',...
    'YMinorTick','on',...
    'LineWidth',1.25)
set(gca, ...
    'FontName','Helvetica');
set(gcf,'color','w','PaperPositionMode','auto')

datetickzoom('x','HH:MM:SS.FFF','keepticks')
linkaxes([ax(1),ax(2),ax(3)],'xy')
%save figure?
prompt = 'Save Figure? [y/n] ';
result = input(prompt,'s');
if strcmp(result,'y') || strcmp(result,'yes');
    fpath = savefigdir;fname = [vname '_ENUtimeseries'];
    export_fig([fpath fname],'-png','-m1','-r900','-opengl')
    disp(['Figure ' fname '.png saved'])
end
if strcmp(result,'n') || strcmp(result,'no');
end