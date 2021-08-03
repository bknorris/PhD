function qp_aqdpdirs(t,u,v,p,intv,fs,lat,lon,name)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% A quick-and-dirty plot of aquadopp velocities using quivers. Velocities
% are depth averaged, then are time averaged based on the 'intv' input.
% Quiver arrows are colored based on the time average interval. This script
% plots two subplots: the Aquadopp's pressure signal as reference for the
% timing of the tide, and the quiver plot.
%
% Inputs:
% t - Aquadopp timebase
% u,v - u and v velocity components
% p - Aquadopp pressure time series
% intv -  interval in MINUTES over which to time-average the data
% fs - the sample rate of the instrument in Hz
% lat,lon - the latitude and longitude of the instrument deployment
% name - a name to label the quiver plot with. Must be a string.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure
f1 = figure(1);
set(f1,'PaperOrientation','landscape',...
    'position',[200 200   1200  600]);
set(gcf,'color','w','PaperPositionMode','auto')    

avt = fs*60*intv;

subplot(121)
%plot pressure record for reference    
idx = [1 avt:avt:length(p)];
c = jet(length(idx));
for j = 1:length(idx)-1
    tx = idx(j):idx(j+1);
    plot(t(tx),p(tx),'Color',c(j,:)), hold on
end
title(['Pressure Signal: ' datestr(t(1),'dd-mm-yyyy')])
datetick('x','HH:MM:SS','keepticks','keeplimits')
hold off
set(gca,'box','on','position',[0.08 0.1 0.38 0.8])
xlabel('Time')
ylabel('Depth (m)')

subplot(122)

%depth average
u = nanmean(u,2);
v = nanmean(v,2);

%time average & plot quivers
for j = 1:length(idx)-1
    t = idx(j):idx(j+1);
    U = nanmean(u(t));
    V = nanmean(v(t));
    hold on
    quiver(lon,lat,U,V,0.001,'Color',c(j,:))
    text(lon,lat,name)
end
set(gca,'box','on','position',[0.54 0.1 0.43 0.8])
xlabel('Longitude (deg)')
ylabel('Latitude (deg)')

