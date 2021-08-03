%mfile to plot barotropic tides on projected map
%need to run get_bar_tide_coef_SoCal first to get the ADCIRC results

load /Projects/PalosVerdes/pvmapping/SoCalBathy
% coast = load('CAcoastline.dat'); - saved in SoCalBathy mat file for quicker loading

% [sites,LAT,LONG] = textread('SoCalSites.txt','%s%f%f');
%get deployment locations
%find appropriate coastline
lonmin = min(min(longs));lonmax = max(max(longs));latmin = min(min(lats));latmax = max(max(lats));
% indx = find(((CAcoastline(:,1)>lonmin)&(CAcoastline(:,1)<lonmax))&((CAcoastline(:,2)>latmin)&(CAcoastline(:,2)<latmax)));
% coast = CAcoastline(indx,:);
indx2 = find(zed>=0);
bathy = zed;bathy(indx2) = NaN;

lonrange = [-118.083 -117.9];latrange = [33.55 33.65];

%initialize the map
m_proj('UTM','zone',11,'hemisphere',0,'long',lonrange,'lat',latrange)
[xMap,yMap]=m_ll2xy(longs,lats);
xM=xMap(1,:); yM=yMap(:,1);
earth = [191/255 191/255 10/255];
figure
hold on
[cs,h] = m_contour(longs,lats,zed,[-1300:100:-100,-60:15:-15]);
set(h,'Color',earth)
h1 = m_line(coast(:,1),coast(:,2));
set(h1,'Color','k');
m_grid
orient landscape


% m_plot(SMB.lon,SMB.lat,'+k')
% m_plot(LA.lon,LA.lat,'+k')
% m_plot(HB.lon,HB.lat,'+k')
% for i = 1:length(LA.sites),m_text(LA.lon(i),LA.lat(i),LA.sites(i)),end
% for i = 1:length(HB.sites),m_text(HB.lon(i),HB.lat(i),HB.sites(i)),end
% for i = 1:length(SMB.sites),m_text(SMB.lon(i),SMB.lat(i),SMB.sites(i)),end


%plot M2 tides (first constituent in list)
tide_vec(10,M651.lon,M651.lat,M651.tide.majU(1),M651.tide.majV(1),M651.tide.minU(1),M651.tide.minV(1),'k')
tide_vec(10,M653.lon,M653.lat,M653.tide.majU(1),M653.tide.majV(1),M653.tide.minU(1),M653.tide.minV(1),'k')
tide_vec(10,M655.lon,M655.lat,M655.tide.majU(1),M655.tide.majV(1),M655.tide.minU(1),M655.tide.minV(1),'k')
tide_vec(10,MHB06.lon,MHB06.lat,MHB06.tide.majU(1),MHB06.tide.majV(1),MHB06.tide.minU(1),MHB06.tide.minV(1),'k')
tide_vec(10,MHB08.lon,MHB08.lat,MHB08.tide.majU(1),MHB08.tide.majV(1),MHB08.tide.minU(1),MHB08.tide.minV(1),'k')
tide_vec(10,MHB10.lon,MHB10.lat,MHB10.tide.majU(1),MHB10.tide.majV(1),MHB10.tide.minU(1),MHB10.tide.minV(1),'k')
tide_vec(10,M658.lon,M658.lat,M658.tide.majU(1),M658.tide.majV(1),M658.tide.minU(1),M658.tide.minV(1),'k')
tide_vec(10,MOCT.lon,MOCT.lat,MOCT.tide.majU(1),MOCT.tide.majV(1),MOCT.tide.minU(1),MOCT.tide.minV(1),'k')
tide_vec(10,MOCQ.lon,MOCQ.lat,MOCQ.tide.majU(1),MOCQ.tide.majV(1),MOCQ.tide.minU(1),MOCQ.tide.minV(1),'k')
tide_vec(10,MOCS.lon,MOCS.lat,MOCS.tide.majU(1),MOCS.tide.majV(1),MOCS.tide.minU(1),MOCS.tide.minV(1),'k')
tide_vec(10,MOCR.lon,MOCR.lat,MOCR.tide.majU(1),MOCR.tide.majV(1),MOCR.tide.minU(1),MOCR.tide.minV(1),'k')

%ADCIRC Results - M2 first
tide_vec(10,loni(1),lati(1),ADC.MAJu(1,4),ADC.MAJv(1,4),ADC.MINu(1,4),ADC.MINv(1,4),'r')
tide_vec(10,loni(2),lati(2),ADC.MAJu(2,4),ADC.MAJv(2,4),ADC.MINu(2,4),ADC.MINv(2,4),'r')
tide_vec(10,loni(3),lati(3),ADC.MAJu(3,4),ADC.MAJv(3,4),ADC.MINu(3,4),ADC.MINv(3,4),'r')
tide_vec(10,loni(4),lati(4),ADC.MAJu(4,4),ADC.MAJv(4,4),ADC.MINu(4,4),ADC.MINv(4,4),'r')
tide_vec(10,loni(5),lati(5),ADC.MAJu(5,4),ADC.MAJv(5,4),ADC.MINu(5,4),ADC.MINv(5,4),'r')
tide_vec(10,loni(6),lati(6),ADC.MAJu(6,4),ADC.MAJv(6,4),ADC.MINu(6,4),ADC.MINv(6,4),'r')
tide_vec(10,loni(7),lati(7),ADC.MAJu(7,4),ADC.MAJv(7,4),ADC.MINu(7,4),ADC.MINv(7,4),'r')
tide_vec(10,loni(8),lati(8),ADC.MAJu(8,4),ADC.MAJv(8,4),ADC.MINu(8,4),ADC.MINv(8,4),'r')
tide_vec(10,loni(9),lati(9),ADC.MAJu(9,4),ADC.MAJv(9,4),ADC.MINu(9,4),ADC.MINv(9,4),'r')
tide_vec(10,loni(10),lati(10),ADC.MAJu(10,4),ADC.MAJv(10,4),ADC.MINu(10,4),ADC.MINv(10,4),'r')
tide_vec(11,loni(11),lati(11),ADC.MAJu(11,4),ADC.MAJv(11,4),ADC.MINu(11,4),ADC.MINv(11,4),'r')

[lhp3,lht3] = m_vec(10,-118.06666,33.63333,5,0,'k','key','Depth Averaged M2 Currents','headlength',0);
[lhp3,lht3] = m_vec(10,-118.06666,33.63,5,0,'r','key','ADCIRC Model','headlength',0);
[lhp3,lht3] = m_vec(10,-118.06666,33.62666,10,0,'k','key','10 CM/S','headlength',0);


clabel(cs,h,'manual','fontsize',10,'FontName','Arial')
