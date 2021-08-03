function  aquadopp=aqdp_read_plot(filename)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Read in currents, backscatter and sen file
%
%  Useage: be sure to run in the same directory as the files to be loaded/
%  plotted. 'filename' is the datafile name without the terminal string.
%
%  Developed by Dr. Julia C Mullarney, University of Waikato, New Zealand 
%  c. 2013
%
%  Additions made by Benjamin K Norris, 2014
%  Edits: yeardate changed to datenum for ease of plotting dates
%         added axis limits on heading plots to represent direction
%         added metadata sub-structure to output file
%         added routine to prompt user for data processing progress plots
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('Loading the data')
[aquadopp.vel_b1]= textread([filename '.v1'],'','delimiter',' /:','headerlines',0);
[aquadopp.vel_b2]= textread([filename '.v2'],'','delimiter',' /:','headerlines',0);
[aquadopp.vel_b3]= textread([filename '.v3'],'','delimiter',' /:','headerlines',0);


[aquadopp.beam1]= textread([filename '.a1'],'','delimiter',' /:','headerlines',0);
[aquadopp.beam2]= textread([filename '.a2'],'','delimiter',' /:','headerlines',0);
[aquadopp.beam3]= textread([filename '.a3'],'','delimiter',' /:','headerlines',0);


[month, day, year, hour, minute, second, aquadopp.error_code, aquadopp.status_code, ...
    aquadopp.batt_voltage, aquadopp.soundspeed, aquadopp.heading, aquadopp.pitch, aquadopp.roll,...
        aquadopp.pressure, aquadopp.temperature, ~, ~]= textread([filename '.sen'],'','delimiter',' /:','headerlines',0);

aquadopp.datenum = datenum(year,month,day,hour,minute,second);
aquadopp.yearday = yearday(year,month,day+hour/24+minute/(24*60)+second/(24*60*60)); %to make Julia happy

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Read in hdr file and check settings/errors
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[hdrfile]= textread([filename '.hdr'],'%s ','delimiter','\t','headerlines',0);

% the following updated for 2014 Fall Mekong Deployment AQDP hdrs
aquadopp.metadata.checksumerrors= str2num(hdrfile{4}(39:end));
aquadopp.metadata.starttime=hdrfile{5}(39:end);
aquadopp.metadata.endtime=hdrfile{6}(39:end);
aquadopp.metadata.dt=str2num(hdrfile{9}(39:40));
aquadopp.metadata.ncells=str2num(hdrfile{10}(39:end));
aquadopp.metadata.cellsize=str2num(hdrfile{11}(39:41))/100;
aquadopp.metadata.measload=str2num(hdrfile{13}(39:41));
aquadopp.metadata.blankdist=str2num(hdrfile{15}(39:42));
aquadopp.metadata.wavemeas=hdrfile{18}(39:end);
aquadopp.metadata.pwrlevel=hdrfile{19}(39:end);
aquadopp.metadata.measint=str2num(hdrfile{20}(39:40));
aquadopp.metadata.spb=str2num(hdrfile{21}(39:end));
aquadopp.metadata.samprate=hdrfile{22}(39:end);
aquadopp.metadata.coordsys=hdrfile{28}(39:end);
aquadopp.metadata.cmt=hdrfile{38}(39:end);
aquadopp.metadata.serialno=str2num(hdrfile{74}(43:end));


if aquadopp.metadata.checksumerrors==0
    disp('There are no errors :)')
end
if aquadopp.metadata.checksumerrors~=0
    disp('There are checksum errors')
end

if strcmp(aquadopp.metadata.wavemeas,'DISABLED')
    disp('Wave measurements are ON')
end

if strcmp(aquadopp.metadata.coordsys,'BEAM')    
    disp('Coordinate system is in BEAM')
end

if strcmp(aquadopp.metadata.coordsys,'XYZ')    
    disp('Coordinate system is in XYZ')
end

if strcmp(aquadopp.metadata.coordsys,'ENU')    
    disp('Coordinate system is in Earth')
end

aquadopp.rangebins= aquadopp.metadata.blankdist+...
    aquadopp.metadata.cellsize:aquadopp.metadata.cellsize:aquadopp.metadata.blankdist+...
    aquadopp.metadata.ncells*aquadopp.metadata.cellsize; %This should now be correct 2014.

prompt = 'Plot figures [y/n]? ';
result = input(prompt,'s');

if strcmp(result,'y') || strcmp(result,'yes');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Basic plots
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure(1)
subplot(211)
    plot(aquadopp.datenum, aquadopp.pressure,'b')
    datetick('x','keepticks')
    title(['Aquadopp ' num2str(aquadopp.metadata.serialno) ': pressure'])
    ylabel('P (dbar)')
subplot(212)
    plot(aquadopp.datenum, aquadopp.temperature,'b')
    datetick('x','keepticks')
    title('temperature')
    ylabel('T (deg C)')
    xlabel('Date')

figure(2)
m1= subplot(311);
    plot(aquadopp.datenum,aquadopp.heading,'r')
    set(m1,'Ylim', [0 360],'YTick',[0:90:360])
    datetick('x','keepticks')
    title(['Aquadopp ' num2str(aquadopp.metadata.serialno) ':heading'])
m2= subplot(312);
    plot(aquadopp.datenum,aquadopp.pitch,'r')
    set(m2,'Ylim', [0 360],'YTick',[0:90:360])
    datetick('x','keepticks')
    title('pitch')
m3 = subplot(313);
    plot(aquadopp.datenum,aquadopp.roll,'r')
    set(m3,'Ylim', [0 360],'YTick',[0:90:360])
    datetick('x','keepticks')
    title('roll')
    xlabel('Date')


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Plot backscatter
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure(3)
ax(1) = subplot(311);
    uimagesc(aquadopp.datenum,aquadopp.rangebins, aquadopp.beam1')
    hold on
    plot(aquadopp.datenum,aquadopp.pressure,'k')
    datetick('x','keepticks')
    title(['Aquadopp ' num2str(aquadopp.metadata.serialno) ': backscatter'])
    c = colorbar;
    axis(c);
    ylabel('Range (m)')
ax(2) = subplot(312);
    uimagesc(aquadopp.datenum,aquadopp.rangebins,aquadopp.beam2')
    hold on
    plot(aquadopp.datenum,aquadopp.pressure,'k')
    datetick('x','keepticks')
    c = colorbar;
    axis(c);
    ylabel('Range (m)')
ax(3) = subplot(313);
    uimagesc(aquadopp.datenum,aquadopp.rangebins, aquadopp.beam3')
    hold on
    plot(aquadopp.datenum,aquadopp.pressure,'k')
    datetick('x','keepticks')
    xlabel('Date')
    c = colorbar;
    axis(c);
    ylabel('Range (m)')
    set(ax,'YDir','Normal','YLim',[0 3])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Plot velocities
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if strcmp(aquadopp.metadata.coordsys,'BEAM') 
    figure(4)
    ax(1) = subplot(311);
        uimagesc(aquadopp.datenum,aquadopp.rangebins, aquadopp.vel_b1')
        hold on
        plot(aquadopp.datenum,aquadopp.pressure,'k')
        datetick('x','keepticks')
        title(['Aquadopp ' num2str(aquadopp.metadata.serialno) ': beam velocities'])
        c = colorbar;
        axis(c);
        ylabel('Range (m)')
    ax(2) = subplot(312);
        uimagesc(aquadopp.datenum,aquadopp.rangebins,aquadopp.vel_b2')
        hold on
        plot(aquadopp.datenum,aquadopp.pressure,'k')
        datetick('x','keepticks')
        c = colorbar;
        axis(c);
        ylabel('Range (m)')
    ax(3) = subplot(313);
        uimagesc(aquadopp.datenum,aquadopp.rangebins, aquadopp.vel_b3')
        hold on
        plot(aquadopp.datenum,aquadopp.pressure,'k')
        datetick('x','keepticks')
        xlabel('Date')
        c = colorbar;
        axis(c);
        ylabel('Range (m)')
        set(ax,'YDir','Normal','YLim',[0 3])
end

if strcmp(result,'n') || strcmp(result,'no');
end
end
   