function  aquadopp=aqdp_read(filename)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Read in currents, backscatter and sen file for Aquadopps
%
%  Useage: be sure to run in the same directory as the files to be loaded/
%  plotted. 'filename' is the datafile name without the terminal string.
%
%  Developed by Dr. Julia C Mullarney, University of Waikato, New Zealand 
%  c. 2013
%
%  Additions made by Benjamin K Norris, 2014
%  Edits: added axis limits on heading plots to represent direction.
%         added metadata sub-structure to output file.
%         added routine to prompt user for data processing progress plots.
%         added routine to extract the instrument transformation matrix
%         from the header file.
%
%  Contains a partial adaptation of READ_aquapro_beam.m 
%  Copyright 2004 
%  USGS Woods Hole Field Center
%  Written by Charlene Sullivan
%  csullivan@usgs.gov
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Read in hdr file and check settings/errors
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp('Loading the data')
[hdrfile]= textread([filename '.hdr'],'%s ','delimiter','\t','headerlines',0);
ind = strfind(hdrfile,'Sampling rate')'; %the header files sometimes have different dimensions
ix = cellfun(@isempty,ind);ind(ix)={nan};
sampind = find([ind{:}]>=1);
samprate = hdrfile{sampind}(39:end); %sample rate

%% Metadata, sensor, velocity and amplitude data are all sample rate specific
if strcmp(samprate,'1 Hz') == 1
    sampmode = 'CONTINUOUS';
    disp('Instrument recorded in CONTINUOUS mode')
else
    sampmode = 'BURST';
    disp('Instrument recorded in BURST mode')
end
switch sampmode
    case 'CONTINUOUS'
        aquadopp.metadata.checksumerrors= str2num(hdrfile{4}(39:end));
        aquadopp.metadata.starttime=hdrfile{5}(39:end);
        aquadopp.metadata.endtime=hdrfile{6}(39:end);
        aquadopp.metadata.dt=str2num(hdrfile{9}(39:40));
        aquadopp.metadata.cellsize=str2num(hdrfile{11}(39:41))/100;
        aquadopp.metadata.measload=str2num(hdrfile{13}(39:41));
        aquadopp.metadata.blankdist=str2num(hdrfile{15}(39:42));
        aquadopp.metadata.wavemeas=hdrfile{18}(39:end);
        aquadopp.metadata.pwrlevel=hdrfile{19}(39:end);
        aquadopp.metadata.measint=str2num(hdrfile{20}(39:40));
        aquadopp.metadata.spb=str2num(hdrfile{21}(39:end));
        aquadopp.metadata.samprate=hdrfile{22}(39:end);
        aquadopp.metadata.coordsys=hdrfile{28}(39:end);
        aquadopp.metadata.salinity=hdrfile{31}(39:42);
        aquadopp.metadata.cmt=hdrfile{38}(39:end);
        aquadopp.metadata.serialno=str2num(hdrfile{74}(43:end));
        ind = strfind(hdrfile,'Transformation matrix')'; %find the transform matrix
        ix = cellfun(@isempty,ind);ind(ix)={nan};
        hdrind = find([ind{:}]>=1);
        t1 = str2num(hdrfile{hdrind}(39:end));
        t2 = str2num(hdrfile{hdrind+1});
        t3 = str2num(hdrfile{hdrind+2});
        aquadopp.transfm = [t1;t2;t3];
        
        %load sensor data
        senData = load([filename '.sen']);
        year = senData(:,3);
        month = senData(:,1);
        day = senData(:,2);
        hour = senData(:,4);
        minute = senData(:,5);
        second = senData(:,6);
        aquadopp.err_code = senData(:,7);    %error code
        aquadopp.sta_code = senData(:,8);    %status code
        aquadopp.batt = senData(:,9);        %Battery (volts)
        aquadopp.sound_spd = senData(:,10);  %soundspeed (m/s)
        aquadopp.heading = senData(:,11);    %Heading (deg)
        aquadopp.pitch = senData(:,12);      %Pitch (deg)
        aquadopp.roll = senData(:,13);       %Roll (deg)
        aquadopp.pressure = senData(:,14);   %Pressure (dBar)
        aquadopp.temperature = senData(:,15);%temperature (degrees C)
        
        % load velocities. these are dimensioned
        % [nBursts x nCells]
        disp('Loading Velocities...')
        aquadopp.vel_b1 = load([filename,'.v1']);disp('Vel1')
        aquadopp.vel_b2 = load([filename,'.v2']);disp('Vel2')
        aquadopp.vel_b3 = load([filename '.v3']);disp('Vel3')
        
        % load amplitudes. these are also dimensioned
        % [nBursts x nCells]
        disp('Loading Amplitudes...')
        aquadopp.beam1 = load([filename '.a1']);disp('Beam1')
        aquadopp.beam2 = load([filename '.a2']);disp('Beam2')
        aquadopp.beam3 = load([filename '.a3']);disp('Beam3')
        [nRows, nCol] = size(aquadopp.vel_b1);
        aquadopp.nBursts = nRows;
        aquadopp.nCells = nCol;
        aquadopp.burst = [1:nRows]';
        
    case 'BURST'
        aquadopp.metadata.checksumerrors= str2num(hdrfile{4}(39:end));
        aquadopp.metadata.starttime=hdrfile{5}(39:end);
        aquadopp.metadata.endtime=hdrfile{6}(39:end);
        aquadopp.metadata.dt=str2num(hdrfile{9}(39:42));
        aquadopp.metadata.cellsize=str2num(hdrfile{10}(39:41))/1000;
        aquadopp.metadata.measload=str2num(hdrfile{22}(39:41));
        aquadopp.metadata.blankdist=str2num(hdrfile{21}(39:42));
        aquadopp.metadata.pwrlevel=hdrfile{30}(39:end);
        aquadopp.metadata.burstsamp=hdrfile{23}(39:end);
        aquadopp.metadata.spb=str2num(hdrfile{24}(39:end));
        aquadopp.metadata.samprate=hdrfile{25}(39:end);
        aquadopp.metadata.coordsys=hdrfile{32}(39:end);
        aquadopp.metadata.salinity=hdrfile{34}(39:42);
        aquadopp.metadata.cmt=hdrfile{41}(39:end);
        aquadopp.metadata.serialno=str2num(hdrfile{79}(43:end));
        ind = strfind(hdrfile,'Transformation matrix')'; %find the transform matrix
        ix = cellfun(@isempty,ind);ind(ix)={nan};
        hdrind = find([ind{:}]>=1);
        t1 = str2num(hdrfile{hdrind}(39:end));
        t2 = str2num(hdrfile{hdrind+1});
        t3 = str2num(hdrfile{hdrind+2});
        aquadopp.transfm = [t1;t2;t3];
        
        %load sensor data
        senData = load([filename '.sen']);
        year = senData(:,3);
        month = senData(:,1);
        day = senData(:,2);
        hour = senData(:,4);
        minute = senData(:,5);
        second = senData(:,6);
        aquadopp.err_code = senData(:,9);    %error code
        aquadopp.sta_code = senData(:,10);   %status code
        aquadopp.batt = senData(:,11);       %Battery (volts)
        aquadopp.sound_spd = senData(:,12);  %soundspeed (m/s)
        aquadopp.heading = senData(:,13);    %Heading (deg)
        aquadopp.pitch = senData(:,14);      %Pitch (deg)
        aquadopp.roll = senData(:,15);       %Roll (deg)
        aquadopp.pressure = senData(:,16);   %Pressure (dBar)
        aquadopp.temperature = senData(:,17);%temperature (degrees C)
        aquadopp.xsen1 = senData(:,18);      %analog input 1
        aquadopp.xsen2 = senData(:,19);      %analog input 2
        
        % load velocities. these are dimensioned [nBursts x nCells]. 
        % Will need to crop out first two columns - sample # and burst
        % counter
        disp('Loading Velocities...')
        aquadopp.vel_b1 = load([filename,'.v1']);aquadopp.vel_b1 = aquadopp.vel_b1(:,3:end);disp('Vel1')
        aquadopp.vel_b2 = load([filename,'.v2']);aquadopp.vel_b2 = aquadopp.vel_b2(:,3:end);disp('Vel2')
        aquadopp.vel_b3 = load([filename '.v3']);aquadopp.vel_b3 = aquadopp.vel_b3(:,3:end);disp('Vel3')

        % load amplitudes. these are also dimensioned [nBursts x nCells].
        % Will need to crop out first two columns
        disp('Loading Amplitudes...')
        aquadopp.beam1 = load([filename '.a1']);aquadopp.beam1 = aquadopp.beam1(:,(3:end));disp('Beam1')
        aquadopp.beam2 = load([filename '.a2']);aquadopp.beam2 = aquadopp.beam2(:,(3:end));disp('Beam2')
        aquadopp.beam3 = load([filename '.a3']);aquadopp.beam3 = aquadopp.beam3(:,(3:end));disp('Beam3')
        
        % load amplitudes. these are also dimensioned [nBursts x nCells].
        % Will need to crop out first two columns
        disp('Loading Correlations...')
        aquadopp.cor1 = load([filename '.c1']);aquadopp.cor1 = aquadopp.cor1(:,(3:end));disp('Cor1')
        aquadopp.cor2 = load([filename '.c2']);aquadopp.cor2 = aquadopp.cor2(:,(3:end));disp('Cor2')
        aquadopp.cor3 = load([filename '.c3']);aquadopp.cor3 = aquadopp.cor3(:,(3:end));disp('Cor3')
        
        [nRows, nCol] = size(aquadopp.vel_b1);
        aquadopp.nBursts = nRows;
        aquadopp.nCells = nCol;
        aquadopp.burst = [1:nRows]';

end

aquadopp.datenum = datenum(year,month,day,hour,minute,second);
aquadopp.yearday = yearday(year,month,day+hour/24+minute/(24*60)+second/(24*60*60)); %to make Julia happy

clear year month day hour minute second senData
if aquadopp.metadata.checksumerrors==0
    disp('There are no errors :)')
end
if aquadopp.metadata.checksumerrors~=0
    disp('There are checksum errors')
end

if strcmp(aquadopp.metadata.coordsys,'BEAM')    
    disp('Coordinate system is in BEAM')
end

if strcmp(aquadopp.metadata.coordsys,'ENU')    
    disp('Coordinate system is in Earth')
end

aquadopp.rangebins= aquadopp.metadata.blankdist+...
    aquadopp.metadata.cellsize:aquadopp.metadata.cellsize:aquadopp.metadata.blankdist+...
    aquadopp.nCells*aquadopp.metadata.cellsize; %This should now be correct 2014.
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
   