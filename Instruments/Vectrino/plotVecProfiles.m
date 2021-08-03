function plotVecProfiles(Data,hcaxis,vcaxis)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%Simple plot function for the Nortek Vectrino Profiler II
% Inputs: 
%       Data: vectrino Data structire
%       hcaxis: scaling for horizontal color axis (0-1)
%       vcaxis = verticalcolor axis specification for velocity color
%        scaling
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%define variables
time = Data.Time;
if isfield(Data,'Profiles_VelBeam1')
    b1 = Data.Profiles_VelBeam1;
    b2 = Data.Profiles_VelBeam2;
    b3 = Data.Profiles_VelBeam3;
    b4 = Data.Profiles_VelBeam4;
    titles = {'Beam 1';'Beam 2';'Beam 3';'Beam 4'};
elseif isfield(Data,'Profiles_VelX')
    b1 = Data.Profiles_VelX;
    b2 = Data.Profiles_VelY;
    b3 = Data.Profiles_VelZ1;
    b4 = Data.Profiles_VelZ2;
    titles = {'X';'Y';'Z1';'Z2'};
end
rangebins = Data.Profiles_Range;
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
ax = tight_subplot(4,1,0.05,[0.1 0.05],[0.1 0.05]);
axes(ax(1))
imagesc(time,rangebins,b1')
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
imagesc(time,rangebins,b2')
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
imagesc(time,rangebins,b3')
title(['\bf\it' titles{3} ' Velocities'])
caxis([min2 max2])
c = colorbar;
cbh = get(c,'Title');
titleString = '[m/s]';
set(cbh ,'String',titleString);
set(c,'YTick',[min2:step2:max2]);
g = gca;
set(g,'XTickLabel',[])
ylabel('\bf\itrange (m)')

axes(ax(4))
imagesc(time,rangebins,b4')
title(['\bf\it' titles{4} ' Velocities'])
caxis([min2 max2])
c = colorbar;
cbh = get(c,'Title');
titleString = '[m/s]';
set(cbh ,'String',titleString);
set(c,'YTick',[min2:step2:max2]);
ylabel('\bf\itrange (m)')
xlabel(['\bf\itTime on ' date])
set(ax,...
    'YTick',[0.04:0.01:0.07],...
    'YLim', [0.04 0.07],...
    'XGrid','on',...
    'GridLineStyle','-',...
    'TickDir','in',...
    'TickLength',[.01 .01],...
    'XMinorTick','off',...
    'YMinorTick','on',...
    'LineWidth',1.25)
set(gca, ...
    'FontName','Helvetica');
set(gcf,'color','w','PaperPositionMode','auto')

datetickzoom('x','HH:MM:SS.FFF','keepticks')
linkaxes([ax(1),ax(2),ax(3),ax(4)],'xy')
end
