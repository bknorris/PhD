clear
% Load and plot the AQUADOPP currents:
% This script takes the Aquadopp files specified in dirc and depth-averages
% the velocity vectors U and V to create quiver plots in order to check the
% orientation of the currents relative to the coastline. It plots a map
% using the Global Shorelines Shapefile as well as outputs .kml files to
% view the data in Google Earth. This script is dependent on a
% 'Instruments.csv' document that contains the Lat/Long data of the
% instrument deployments under the headings:
% Site/Instrument/Lat/Long

% Note: Before using this script, make sure to use
% aqdp_dataprocess.m to process the raw data, and burstaverage.m
% to generate the time-averaged velocities.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% EDIT THESE VALUES
date = '30092014';
dirc = 'c:\Users\bkn5\Projects\Mekong_F2014\Data\Aquadopp\F2F\';
sca = 0.1; %map scalar

yyyy = 2014; mm = 09; dd = 29;hr = 07; mi = 00; sec = 00;
start = yearday(yyyy,mm,dd+hr/24+mi/(24*60)+sec/(24*60*60));

dd = 29;hr = 19; mi = 00; sec = 00;
finish = yearday(yyyy,mm,dd+hr/24+mi/(24*60)+sec/(24*60*60));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
sdir = 'c:\Users\bkn5\Projects\Mekong_F2014\ArcGis\Landboundaries\GSHHS_shp\f\';
dirll = 'c:\Users\bkn5\Projects\Mekong_F2014\Documents\Deployments\';
dirf = ['c:\Users\bkn5\Projects\Mekong_F2014\Figures\AQDP\' date '\'];
fname = 'Instruments.csv';
l = m_shaperead([sdir 'GSHHS_f_L1_cld']);
files = dir([dirc '*_p.mat']);
filelist = {files.name};

for ii = 1:length(filelist)
    varname = filelist{ii};
    newvarname = regexprep(varname,'_f.mat','');
    load([dirc varname])
    
    yearday = aqdp.yearday;
    rangebins = aqdp.rangebins;
    spb = aqdp.metadata.spb;
    nBursts = length(aqdp.burst);
    hrs = (yearday-(floor(yearday(1))))*24;
    burstt = 1:nBursts;
    hrs = hrs(burstt);
    yearday = yearday(burstt);
    dindx = find(yearday >= start & yearday <= finish);
    pressure = aqdp.pressure(burstt);pressure = pressure(dindx);
    names = {'MA-1';'MA-2a';'MA-2b';'Control';'MA-3';'Flats2Forest';'MA-4'};
    if strcmp(date,'24092014') == 1;
        siteno = 1;
    end
    if strcmp(date,'25092014') == 1;
        siteno = 2;
    end
    if strcmp(date,'26092014') == 1;
        siteno = 3;
    end
    if strcmp(date,'27092014') == 1;
        siteno = 4;
    end
    if strcmp(date,'28092014') == 1;
        siteno = 5;
    end
    if strcmp(date,'29092014') == 1;
        siteno = 6;
    end
    if strcmp(date,'30092014') == 1;
        siteno = 6;
    end
    if strcmp(date,'02102014') == 1;
        siteno = 7;
    end
    if strcmp(date,'03102014') == 1;
        siteno = 7;
    end
    if strcmp(newvarname,'HR5116') == 1;
        newvarname = regexprep(newvarname,'HR','');
    end
    if strcmp(newvarname,'HR5117') == 1;
        newvarname = regexprep(newvarname,'HR','');
    end
    
    %Get instruments Lat/Long from the Locations file
    fid = fopen([dirll fname]);
    str = fgetl(fid);
    header = textscan(str,'%s','Delimiter',',');
    fields = header{1};
    FRMT = '%s%s%s%s';
    thedata = textscan(fid,FRMT,'delimiter',',');
    ind = strfind(thedata{1},names{siteno})';
    ix = cellfun(@isempty,ind);ind(ix)={nan};
    ind2 = find([ind{:}]>=1);
    ind3 = strfind(thedata{2}(ind2),newvarname);
    ix = cellfun(@isempty,ind3);ind3(ix)={nan};
    ind4 = find([ind3{:}]>=1);
    
    %Compute depth averages of burst-averaged data
    [a,b] = size(aqdp.uav);%     lat = cell2mat(thedata{3}(ind2(ind4)));lon = cell2mat(thedata{4}(ind2(ind4)));
%     lat = regexp(lat, '(\d+?)*\d+(\.\d*)?', 'match');
%     lon = regexp(lon, '(\d+?)*\d+(\.\d*)?', 'match');
%     dd = str2num(lat{1});mm = str2num(lat{2}); ss = str2num(lat{3});
%     lat = [dd mm ss];
%     dd = str2num(lon{1});mm = str2num(lon{2}); ss = str2num(lon{3});
%     lon = [dd mm ss];

    for i = 1:a
        ud(i,:) = nanmean(aqdp.uav(i,:));
        vd(i,:) = nanmean(aqdp.vav(i,:));
    end
%     ud = ud(dindx);
%     vd = vd(dindx);
    [a,b] = size(ud);
    lat1 = [9.491309,9.501882,9.512472,9.522294];
    lon1 = [106.24435,106.244113,106.243944,106.2441];
    lat2 = lat1(ii);
    lon2 = lon1(ii);
    latllim = lat2-sca;latulim = lat2+sca;
    lonllim = lon2-sca;lonulim = lon2+sca;
    goods = find(~isnan(ud));
    ud = ud(goods);vd = vd(goods);hrs = hrs(goods);pressure = pressure(goods);
    lat = repmat(lat2,length(goods),b);
    lon = repmat(lon2,length(goods),b);
    c = {'ffff1100';'ffff9900';'ff0099ff';'ff0011ff'}; %create new color for each iteration
    
    % Comparison Map Plot
    f1 = figure(1);
    set(f1,'PaperOrientation','portrait',...
        'position',[100 100   600   600]);
    for k=1:length(l.ncst),
        line(l.ncst{k}(:,1),l.ncst{k}(:,2),'Color','k','LineWidth',1.5);
    end;
    hold on
    axis([lonllim lonulim latllim latulim])
    h(1) = quiver(lon,lat,ud,vd,0.1);hold on
    d(1) = line(lon2,lat2,'marker','square','markersize',5,'color','r');
    d(2) = text(lon2+(sca/10),lat2,names{siteno},'vertical','bottom');
    set(d(2),'FontWeight','bold')
    axis square
    g = gca;
    set(g,...
        'XGrid','on',...
        'YGrid','on',...
        'YMinorTick','on',...
        'XMinorTick','on',...
        'GridLineStyle','--',...
        'YGrid','on',...
        'TickDir','in',...
        'XColor',[.3 .3 .3], ...
        'YColor',[.3 .3 .3], ...
        'LineWidth',1,...
        'TickLength',[.02 .02],...
        'XTick',lonllim:(sca/2):lonulim,...
        'YTick',latllim:(sca/2):latulim)
    box on
    title([varname ' ' date])
    set(gcf,'color','w','PaperPositionMode','auto');
    
    % Compass Plot
    % f1 = figure(2);
    % set(f1,'PaperOrientation','portrait',...
    %     'position',[200 200   600   600]);
    % m1 = compass(ud,vd);
    % title('U,V velocity components - 90 is due North')
    
    
%     % Vector Plots
%     f1 = figure(3);
%     set(f1,'PaperOrientation','portrait',...
%         'position',[200 200   600   600]);
%     subplot(211)
%     xx = 1:length(hrs);
%     yy = zeros(size(xx));
%     line(xx,yy,'LineWidth',1.5,'Color','k')
%     hold on
%     vecplot2(hrs,ud,vd,[hrs(1) hrs(end)],3);
%     ylabel('\bf\it(m/s)')
%     axis([hrs(1) hrs(end) -0.15 0.15])
%     t1 = title([newvarname ' ' date]);
%     set(t1,'FontWeight','bold','FontAngle','italic')
%     g = gca;
%     set(g,'XTick',floor(hrs(1)):1:floor(hrs(end)),...
%         'YTick',-0.15:0.05:0.15);
%     box on
%     subplot(212)
%     g1 = plot(hrs,pressure);
%     set(g1,'LineWidth',2);
%     xlabel('\bf\itHours After the Start of the Experiment')
%     ylabel('\bf\itdbar')
%     axis([hrs(1) hrs(end) 0 1.5])
%     g = gca;
%     set(g,'XTick',floor(hrs(1)):1:floor(hrs(end)));
%     box on
%     set(gcf,'color','w','PaperPositionMode','auto');
%     export_fig([dirf newvarname '_vec'],'-png','-m1','-r900','-opengl')
% 
%     %Create some cool layers for Google Earth
%     kmlstr.f1 = ge_quiver(lon,lat,ud,vd,...
%         'lineColor',c{ii},...
%         'lineWidth',1.2,...
%         'altitudeMode','clampToGround',...
%         'magnitudeScale',0.001,...
%         'msgToScreen',true);
%     kmlstr.f2 = ge_point(lon2,lat2,1,...
%         'altitudeMode','clampToGround',...
%         'msgToScreen',true);
%     kmlstr.f3 = ge_text(lon2,lat2,1,[newvarname '_' date],...
%         'altitudeMode','clampToGround',...
%         'msgToScreen',true);
%     files = {[newvarname '_quiver.kml'];[newvarname '_point.kml'];[newvarname '_text.kml']};
%     f = fieldnames(kmlstr);
%     for i = 1:length(f)
%         ge_output([dirf files{i}],kmlstr.(f{i}));
%     end
%     pause(4)
%     close all
end
disp('Finished!')
clear all



