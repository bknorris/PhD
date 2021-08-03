function ADV = readADVdata(filename,metadata)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%A function to read in raw Nortek Vector single point current meter data
%and create a structure (adv) which contains the metadata, .sen file data,
%and the .dat file data. NOTE: This script will only work for Vectors that
%recorded continuously, as instruments that burst-sampled have their data
%files arranged differently.
%
%   filename - the base filename without the terminal string
%   metadata - a structure of metadata, user defined. See adv_dataprocess.m
%   for an example of the required inputs
%
%   Dependencies: readVEChdr.m
%
%   Contains an adaptation of vec2cdfBurstNative.m
%   Written by Kurt J Rosenberger
%   krosenberger@usgs.gov
%   USGS Pacific Coastal Marine Science Center
%
%   This script was compiled by Benjamin K Norris, 2015
%   University of Waikato, New Zealand
%
%   Updates:
%   11052015: Added warning messages to be displayed if the instrument status
%             code returns out of range values for pitch/roll. Added
%             adjustment to velocities if status reports scaling.
%  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if nargin<1
    help(mfilename)
    return
end
thefillvalue = NaN; %define a fill value for non-elements
instmeta = readVEChdr(filename);

if ~strcmp(instmeta.burstint,'CONTINUOUS')
    %if the instrument did not record in continuous mode
    disp('Error - instrument recorded in BURST mode')
    return
end
%define start/stop times
deploy = datenum(metadata.deployment_date);
recover = datenum(metadata.recovery_date);

%load the data, starting with the sensor data to get time stamps. We cannot
%assume that every burst has the right amount of data, so we'll read the
%.dat file line by line
senFile = strcat(filename,'.sen');
if ~exist(senFile,'file')
    disp('Error - no .sen file is found in this directory')
    return
end

disp('Loading external sensor data')
sendata = load(senFile);
%Find breaks in the timebase to delineate bursts. The sensor data is
%collected at 1Hz, regardless of the velocity sampling rate.
yy = instmeta.senflds.Year;mo = instmeta.senflds.Month;
dd = instmeta.senflds.Day;hh = instmeta.senflds.Hour;
mi = instmeta.senflds.Minute;ss = instmeta.senflds.Second;
dn = datenum([sendata(:,yy) sendata(:,mo) sendata(:,dd) sendata(:,hh) sendata(:,mi) sendata(:,ss)]);

%gooddata should be used to crop datafiles to start and stop times as
%defined as metadata.deployment_date and metadata.recovery_date.
gooddata = find(dn>=deploy&dn<=recover); %indices of good data
datadn = dn(gooddata);

%plot the heading to check that the given start and stop times are actually
%correct
heading = instmeta.senflds.Heading;heading = sendata(:,heading);
heading = heading(gooddata);
disp('USER check that the supplied deployment and recovery times are correct')
f1 = figure(1);
set(f1,'PaperOrientation','portrait',...
    'position',[400 300   800   550]);
plot(datadn,heading) 
ylabel('Degrees')
datetick('x','dd HH:MM:SS','keepticks')
title('Instrument Heading')
set(gca,'XGrid','on','Xlim',[min(datadn) max(datadn)]);
hold on
pause(5)
prompt = 'Are these times correct [y/n]? ';
result = input(prompt,'s');
if strcmp(result,'y') || strcmp(result,'yes');
elseif strcmp(result,'n') || strcmp(result,'no');
    figure(1)
    disp('USER select start/end points where the heading was consistent')
    clear q
    q = ginput(2);
    ind1 = q(1);ind2 = q(2); %new start/end times
    gooddata = find(dn>=ind1&dn<=ind2);
    datadn = dn(gooddata);
end

close(f1)
nsamp = length(gooddata);
startind = gooddata(1);
endind = gooddata(end);
disp(['There are ',num2str(nsamp),' good samples in the file'])

%Crop sensor data to gooddata, load into a structure
senfields = fieldnames(instmeta.senflds);
SENdata = sendata(startind:endind,:);
for i = 8:16
    senfld = senfields{i};
    SEN.(senfld) = SENdata(:,i);
end
SEN.Datetime = datadn;
N = length(SEN.Datetime);fields = fieldnames(SEN);
SEN.int = instmeta.burstint;

%Check status code from each sample to determine instrument
%parameters/potential errors
orient = zeros(N,1);scaling = zeros(N,1);
pitcherr = zeros(N,1);rollerr = zeros(N,1);

for i = 1:N;
    status = num2str(SEN.Statuscode(i));
    orient(i,:) = str2double(status(7));
    scaling(i,:) = str2double(status(6));
    pitcherr(i,:) = str2double(status(5));
    rollerr(i,:) = str2double(status(4));
end

if any(orient)
    disp('Instrument claims it points UP')
else
    disp('Instrument claims it points DOWN')
end
if any(scaling)
    disp('Scaling is set to 0.1mm/s')
else
    disp('Scaling is set to mm/s')
end
if any(pitcherr)
    disp('WARNING: Pitch measurements are out of range')
end
if any(rollerr)
    disp('WARNING: Roll measurements are out of range')
end
clearvars sendata

% put in fill values    
for i = 1:length(SEN)
    fld = fields{i};
    if length(SEN.(fld) == N); %#ok<*ISMT> %i.e. only work on timeseries
        if strcmp(fld,{'time';'time2'}),
            continue
        else
            nanind = isnan(SEN.(fld));SEN.(fld)(nanind) = thefillvalue;
        end
    end
end

%define metadata fields
metadata.start_time = datestr(SEN.Datetime(1),'yyyy-mm-dd HH:MM:SS:FFF');
metadata.stop_time = datestr(SEN.Datetime(end),'yyyy-mm-dd HH:MM:SS:FFF');
metadata.sampleinds = gooddata;
metadata.inst_type = 'Nortek Vector Acoustic Doppler Velocimeter';
metadata.instmeta = instmeta;
metadata.nsamp = nsamp;

%create new burst timebase
%if the first timestamp in the .SEN file = the first timestamp
%of datadn, it is the first sample and the new time series needs to start on
%the second sample in datadn. Continuous datasets can be considered
%as a single burst, so if the above relationship is not true, the timestamp
%needs to start on the first value of datadn.
samprate = instmeta.samprate;
dTs = 1/samprate;
dTsdn = datenum(0,0,0,0,0,dTs); %convert seconds/sample to a serial datenumber
SEN.time = floor(datadn);
SEN.time2 = (datadn - floor(datadn));
time = zeros(samprate,nsamp);
time2 = zeros(samprate,nsamp);
if dn(1) == datadn(1)
     samptime = 0:dTsdn:dTsdn*(samprate-1);
     samptime = samptime+datenum(0,0,0,0,0,1); %add one second of time offset for the start of the burst
     for i = 1:nsamp
         time(:,i) = repmat(SEN.time(i),samprate,1);
         time2(:,i) = SEN.time2(i) + samptime;
     end
     time = reshape(time,nsamp*samprate,1); %create date vector
     time2 = reshape(time2,nsamp*samprate,1); %create time vector
elseif dn(1) ~= datadn(1)
    samptime = 0:dTsdn:dTsdn*(samprate-1);
    for i = 1:nsamp
        time(:,i) = repmat(SEN.time(i),samprate,1);
        time2(:,i) = SEN.time2(i) + samptime;
    end
    time = reshape(time,nsamp*samprate,1); %create date vector 
    time2 = reshape(time2,nsamp*samprate,1); %create time vector
end
dattime = time+time2;    

%now get the data from the .dat file and write incrementally
n = nsamp*samprate;
datid = fopen([filename '.dat']);
count = 1;
if startind == 1 %if startind = 1, begin reading file at first line
elseif startind>1
    while count<((startind*samprate)-(samprate-1)) 
        str = fgetl(datid);count = count+1; %#ok<NASGU> %set the file position counter to the value of count
    end
end

%we can't assume that each burst has the right number of datapoints, so
%we'll set up a structure that does have the right number, and fill them in
%line by line, and any missing samples will get a fill values.
%fields in the output structure:
DAT = struct('Burst',NaN(n,1),'Ens',NaN(n,1),'V1',NaN(n,1),'V2',NaN(n,1),...
    'V3',NaN(n,1),'Amp1',NaN(n,1),'Amp2',NaN(n,1),'Amp3',NaN(n,1),...
    'Snr1',NaN(n,1),'Snr2',NaN(n,1),'Snr3',NaN(n,1),'Cor1',NaN(n,1),...
    'Cor2',NaN(n,1),'Cor3',NaN(n,1),'Pres',NaN(n,1),'Analog1',NaN(n,1),...
    'Analog2',NaN(n,1),'Chk',NaN(n,1));
flds = fieldnames(DAT);
tic
disp('Now reading from .dat file')
for i = 1:n
    str = fgetl(datid);buf = textscan(str,'%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f');
    for j = 1:length(flds)
        fld = flds{j};
        DAT.(fld)(i,:) = buf{j};
    end
    if (i/10000)-floor(i/10000)==0 %count samples read
        disp([num2str(i) ' samples of ' num2str(n) ' read']),
    end
    if rem(i,n)==0
        disp(['Finished reading in ',num2str(i),' samples'])
    end
end

%assume velocities are usually measured in m/s, if scaling is set to
%0.1mm/s, scale to mm/s
if any(scaling)
    DAT.VEL1 = DAT.V1.*10;DAT.VEL2 = DAT.V2.*10;DAT.VEL3 = DAT.V3.*10;
end

%define time vectors for the DAT structure
DAT.datetime = dattime;
dvec = datevec(dattime); %create yearday vector for Julia
yy = dvec(:,1);mo = dvec(:,2);dd = dvec(:,3);
hh = dvec(:,4);mi = dvec(:,5);ss = dvec(:,6);
DAT.yearday = yearday(yy,mo,dd+hh/24+mi/(24*60)+ss/(24*60*60));
DAT.time = time;
DAT.time2 = time2;

% now do the fill values for the RAW data
N = length(DAT.time);flds = fieldnames(DAT);
for i = 1:length(DAT)
    fld = flds{i};
    if length(DAT.(fld) == N); %i.e. only work on timeseries
        if strcmp(fld,{'time';'time2'}),
            continue
        else
            nanind = isnan(DAT.(fld));DAT.(fld)(nanind) = thefillvalue;
        end
    end
end

%copy files to the structure ADV
flds = fieldnames(SEN); %ignore time, time2 and sampint
for i = 1:10
    ADV.Sensor.(flds{i}) = SEN.(flds{i});
end
flds = fieldnames(DAT);
for i = 1:length(flds)
    ADV.(flds{i}) = DAT.(flds{i});
end

disp(['Finished reading and writing ',num2str(n),' samples in ',num2str(toc/60),' minutes'])

%check to make sure we have the height of the pressure sensor. The distance
%of the sensor from the probe is constant for fixed head Vectors. If the
%vector type has a flexible head, we cannot determine this value as easily.
if ~isfield(metadata,'probe_height')
    if strcmp(metadata.orientation,'HORIZONTAL')
        if strcmp(metadata.vector_type,'FIXED')
            metadata.probe_height = metadata.pressure_sensor_height;
        else
            disp('No probe height found in the metadata')
            metadata.probe_height = input('Enter height of the Vector probe, center transducer: ');
        end
    end
end
    
if ~isfield(metadata,'pressure_sensor_height')
    if strcmp(metadata.vector_type,'FIXED')
        switch metadata.orientation
            case 'UP'           %the probe is facing up, pressure sensor below it
                metadata.pressure_sensor_height = metadata.probe_height - 214;
                disp('No pressure sensor height in the metadata, thus it is assumed to be a fixed vector')
                disp('---The orientation is UP, thus the pressure sensor height is set to the probe height - 214mm')
            case 'DOWN'        %the probe is facing down, pressure sensor above it
                metadata.pressure_sensor_height = metadata.probe_height + 214;
                disp('No pressure sensor height in the metadata, thus it is assumed to be a fixed vector')
                disp('---The orientation is DOWN, thus the pressure sensor height is set to the probe height + 214mm')
            case 'HORIZONTAL'
                metadata.pressure_sensor_height = metadata.probe_height;
                disp('No pressure sensor height in the metadata, thus it is assumed to be a fixed vector')
                disp('---The orientation is HORIZONTAL, thus the pressure sensor height is set to the probe height')
        end
    end
    if strcmp(metadata.vector_type,'FLEXIBLE')
        disp('No pressure sensor height in the metadata')
        metadata.pressure_sensor_height = input('Enter height of the Vector case: ');
    end
end
ADV.Metadata = metadata;
clearvars SEN DAT time time2 yy mo dd hh mi ss dattime metadata instmeta orient status pitcherr rollerr
end

