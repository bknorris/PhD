function [SEN,RAW] = vec2cdfBurstNative(baseFile,metadata)
%function aqd2cdf('baseFile','cdfFile','metaFile')
%
%   Function to read in raw Nortek Vector Single point current meter 
%   data and create  a cdf file of the raw data. 
%
%   baseFile - base filename for raw Nortek Vector data (.dat)
%
%   metadata - structure of metadata. See runVec2cdfNativeExample.m for
%   required inputs
%
%       

if nargin<1
    help(mfilename)
    return
end
metadata.history = [];
doublefill = 1e35;
shortfill = -32768;

% get the metadata about setup and log from deployment
files = dir;
for i = 1:length(files)
    names{i} = files(i).name;
end
if ~isempty(strmatch([baseFile '.hdr'],names))
    instmeta = readVECHeader(baseFile);
else
end
deploy = datenum2julian(datenum(metadata.Deployment_date));
recover = datenum2julian(datenum(metadata.Recovery_date));
%now load data  - first load the external sensor data to get time stamps,
%then read in the data from the .DAT file. We can't assume that every burst
%has the right amount of data, so we'll have to go line by line...
senFile = strcat(baseFile,'.sen');
if ~exist(senFile,'file')
    disp('Error - no .sen file found in this directory')
    return
end

disp('Loading external Sensor Data first')
sendata = load(senFile);[n,m] = size(sendata);
%find breaks in timebase to delineate bursts. The sensor
%data is collected at 1Hz, regardless of velocity sampling rate.
jd = julian([sendata(:,3) sendata(:,1) sendata(:,2) sendata(:,4) sendata(:,5) sendata(:,6)]);
DJD = find(diff(jd)>(10/86400))+1;                  %the start indeces of all the bursts
BSTJD = jd(DJD);                                    %the start times of all the bursts
gooddata = find(jd>=deploy&jd<=recover);            %the indeces of all the good data
goodbursts = find(BSTJD>=deploy&BSTJD<=recover);    %the good burst numbers
burstjd = BSTJD(goodbursts);                        %the start times of the good bursts
startind = DJD(goodbursts);                         %the start indeces of the good bursts
%check to see if it was pulled out of the water mid-burst, very likely
spb = instmeta.VECSamplesPerBurst;
if (gooddata(end)-startind(end))<spb
    goodbursts = goodbursts(1:end-1);burstjd = burstjd(1:end-1);startind = startind(1:end-1);
end
if isfield(metadata,'good_bursts')
    disp(['User instructed to process bursts ',num2str(metadata.good_bursts(1)),' to ',num2str(metadata.good_bursts(2))])
    ind1 = near(goodbursts,metadata.good_bursts(1),1);
    ind2 = near(goodbursts,metadata.good_bursts(2),1);
    burstjd = burstjd(ind1:ind2);
    startind = startind(ind1:ind2);
    goodbursts = metadata.good_bursts(1):1:metadata.good_bursts(2);
    nbursts = length(goodbursts);
else
    nbursts = length(goodbursts);
    disp(['There are ',num2str(nbursts),' good bursts in the file'])
end
    
% need to grab the error and status codes from the first good burst for 
% orientation and scaling information
senid = fopen(senFile);
for i = 1:startind(1)+1
    str = fgetl(senid);
end
buf = textscan(str,'%s%s%s%s%s%s%s%s%f%f%f%f%f%f%f%f%f');
error = char(buf{7});status = char(buf{8});
orient = (status(8));scaling = (status(7));
fclose(senid);

senflds = {'Month';'Day';'Year';'Hour';'Minute';'Second';'Error_code';'Status_code';...
    'battery';'soundspeed';'heading';'pitch';'roll';'temperature';'Analog_input';'Checksum'};
    
%now want to average the sensor data by burst to get average pitch, roll
%and heading for transformations
% TODO - save the actual burst ancillary sensor data in a separate file
% with it's own time base, since it is only one Hz
SENDATA = zeros(nbursts,m);
for i = 1:nbursts-1
    p = startind(i);q = startind(i+1);
    SENDATA(i,:) = mean(sendata(p:q,:));
end
enddata = DJD(goodbursts(end)+1);
SENDATA(end,:) = mean(sendata(startind(end):enddata,:));
% SENDATA(end,:) = mean(sendata(startind(end):startind(end)+spb,:));
for i = 9:15
    senfld = senflds{i};
    SEN.(senfld) = SENDATA(:,i);
end
%perform time correction if specified
if isfield(metadata,'TimeCorrection')
    newjd = burstjd + metadata.TimeCorrection/24;
    SEN.jd = newjd;
    metadata.history = [metadata.history 'Time has been corrected by ',num2str(metadata.TimeCorrection),' hours'];
else
    SEN.jd = burstjd;
end
N = length(SEN.jd);flds = fieldnames(SEN);

SEN.int = instmeta.VECBurstInterval/86400;
%grab the pitch roll and heading from another SEN file if requested
if isfield(metadata,'external_compass')
    disp('Looking for sensor data from external instrument for transformations')
    EXTdata = load(metadata.external_compass);
    [pth,nm,ext] = fileparts(metadata.external_compass);
    hdrfile = [pth '/' nm];
    [EXTmeta,SN] = readNortekGenericHeader(hdrfile);
    
    if isfield(EXTmeta,'ProfileInterval')      %i.e., this is a profiler in profile mode, not burst
        EXTjd = julian([EXTdata(:,3) EXTdata(:,1) EXTdata(:,2) EXTdata(:,4) EXTdata(:,5) EXTdata(:,6)]);
        EXTint = EXTmeta.ProfileInterval/86400;
        overlap = find(EXTjd>=SEN.jd(1)&EXTjd<=SEN.jd(end));
        EXT = EXTdata(overlap,:);
        EXTjd = julian([EXT(:,3) EXT(:,1) EXT(:,2) EXT(:,4) EXT(:,5) EXT(:,6)]);
        if length(EXTjd)~=N
            difflen = abs(N-length(EXTjd));
            EXTJD = zeros(N,1);EXTpitch = zeros(N,1);EXTheading = zeros(N,1);EXTroll = zeros(N,1);
            if EXTjd(1)>SEN.jd(1)
                EXTJD(1:difflen) = SEN.jd(1:difflen);EXTJD(difflen+1,end) = EXTjd(1:end);
                EXTheading(1:difflen) = SEN.heading(1:difflen);EXTheading(difflen+1:end) = EXT(SN.Heading,1:end);
                EXTpitch(1:difflen) = SEN.pitch(1:difflen);EXTpitch(difflen+1:end) = EXT(SN.Pitch,1:end);
                EXTroll(1:difflen) = SEN.roll(1:difflen);EXTroll(difflen+1:end) = EXT(SN.Roll,1:end);
            elseif SEN.jd(end)>EXTjd(end)
                EXTJD(difflen:end) = SEN.jd(difflen:end);EXTJD(difflen+1,end) = EXTjd(1:end);
                EXTheading(difflen:end) = SEN.heading(1:difflen);EXTheading(1:difflen-1) = EXT(SN.Heading:end);
                EXTpitch(difflen:end) = SEN.pitch(difflen:end);EXTpitch(1:difflen-1) = EXT(SN.Pitch:end);
                EXTroll(difflen:end) = SEN.roll(difflen:end);EXTroll(1:difflen-1) = EXT(SN.Roll:end);
            end
        end
        if SEN.int~=EXTint
            EXT.heading = interp1(EXTheading,SEN.jd);EXT.roll = interp1(EXTroll,SEN.jd);
            EXT.pitch = interp1(EXTpitch,SEN.jd);EXT.jd = SEN.jd;
        else
            EXT.jd = EXTJD;EXT.heading = EXTheading;EXT.pitch = EXTpitch;EXT.roll = EXTroll;
        end
    elseif isfield(EXTmeta,'BurstInterval')      %i.e., this is a profiler in burst mode
        EXTint = EXTmeta.BurstInterval/86400;
        jddummy = julian([EXTdata(:,3) EXTdata(:,1) EXTdata(:,2) EXTdata(:,4) EXTdata(:,5) EXTdata(:,6)]);
        overlap = find(jddummy>=SEN.jd(1)&jddummy<=SEN.jd(end));
        EXTcrop = EXTdata(overlap,:);
%         EXTjd = julian([EXTcrop(:,3) EXTcrop(:,1) EXTcrop(:,2) EXTcrop(:,4) EXTcrop(:,5) EXTcrop(:,6)]);
%         EXTDJD = find(diff(EXTjd)>(10/86400))+1;                  %the start indeces of all the bursts
%         EXTBSTJD = jd(EXTDJD);                                    %the start times of all the bursts
        [l,o] = size(EXTcrop);nextbst = floor(l/EXTmeta.SamplesPerBurst);ntot = nextbst*EXTmeta.SamplesPerBurst;
        EXTCROP = EXTcrop(1:ntot,:);
%         ext = reshape(EXTcrop,EXTmeta.SamplesPerBurst,nextbst,o);
        extrshp = squeeze(mean(reshape(EXTCROP,EXTmeta.SamplesPerBurst,nextbst,o)));
        EXTjd = julian([extrshp(:,3) extrshp(:,1) extrshp(:,2) extrshp(:,4) extrshp(:,5) extrshp(:,6)]);
        if length(EXTjd)~=N
            difflen = abs(N-length(EXTjd));x = 1:1:N;
            if length(EXTjd)>N
                EXT.heading = interp1(extrshp(:,SN.Heading),x);EXT.pitch = interp1(extrshp(:,SN.Pitch),x);
                EXT.roll = interp1(extrshp(:,SN.Roll),x);EXT.jd = SEN.jd;
            elseif N>length(EXTjd)
                EXTJD = zeros(N,1);EXTpitch = zeros(N,1);EXTheading = zeros(N,1);EXTroll = zeros(N,1);
                if EXTjd(1)>SEN.jd(1)
                    EXTJD(1:difflen) = SEN.jd(1:difflen);EXTJD(difflen+1,end) = EXTjd(1:end);
                    EXTheading(1:difflen) = SEN.heading(1:difflen);EXTheading(difflen+1:end) = extrshp(SN.Heading,1:end);
                    EXTpitch(1:difflen) = SEN.pitch(1:difflen);EXTpitch(difflen+1:end) = extrshp(SN.Pitch,1:end);
                    EXTroll(1:difflen) = SEN.roll(1:difflen);EXTroll(difflen+1:end) = extrshp(SN.Roll,1:end);
                elseif SEN.jd(end)>EXTjd(end)
                    EXTJD(difflen:end) = SEN.jd(difflen:end);EXTJD(difflen+1,end) = EXTjd(1:end);
                    EXTheading(difflen:end) = SEN.heading(1:difflen);EXTheading(1:difflen-1) = extrshp(SN.Heading:end);
                    EXTpitch(difflen:end) = SEN.pitch(difflen:end);EXTpitch(1:difflen-1) = extrshp(SN.Pitch:end);
                    EXTroll(difflen:end) = SEN.roll(difflen:end);EXTroll(1:difflen-1) = extrshp(SN.Roll:end);
                end
                if SEN.int~=EXTint
                    EXT.heading = interp1(EXTheading,x);EXT.roll = interp1(EXTroll,x);
                    EXT.pitch = interp1(EXTpitch,x);EXT.jd = SEN.jd;
                else
                    EXT.jd = EXTJD;EXT.heading = EXTheading;EXT.pitch = EXTpitch;EXT.roll = EXTroll;
                end
            end
        end
    end
    SEN.heading = EXT.heading;SEN.pitch = EXT.pitch;SEN.roll = EXT.roll;
end

        
% put in fill values    
for i = 1:length(SEN)
    fld = flds{i};
    if length(SEN.(fld) == N); %i.e. only work on timeseries
        if strcmp(fld,{'time';'time2'}),
            continue
        else
            nanind = isnan(SEN.(fld));SEN.(fld)(nanind) = doublefill;
        end
    end
end


SEN.time = floor(burstjd); SEN.time2 = (burstjd - floor(burstjd))*86400000;
metadata.goodbursts = goodbursts;

% be very clear about magvar and if it was applied
if ~isfield(metadata,'magnetic_variation_at_site'),
    metadata.magnetic_variation_at_site = metadata.magnetic_variation;
else
    metadata.magnetic_variation_at_site = 0;
end
if isfield(metadata,'magnetic_variation'),
    metadata = rmfield(metadata,'magnetic_variation');
end
metadata.INST_TYPE = 'Nortek Vector acoustic Doppler velocimeter';
% SEN.metadata = metadata;
metadata.instmeta = instmeta;
metadata.nbursts = nbursts;
metadata.burstlength = spb;
metadata.start_time = gregorian(SEN.jd(1));
metadata.stop_time = gregorian(SEN.jd(end));

%create new burst timebase
time = zeros(spb,nbursts);time2 = zeros(spb,nbursts);
dTb = 1/instmeta.VECSamplingRate*1000;
bsttime = 0:dTb:dTb*(spb-1);
for i = 1:nbursts
    time(:,i) = repmat(SEN.time(i),1,spb);
    time2(:,i) = SEN.time2(i) + bsttime;
end

%% define the file
if ~isfield(metadata,'filename')
    cdfFilename = strcat([metadata.instrument_number 'vec-b.cdf']);
%     metadata.ncfilename	= cdfFilename;				% output name including search path
else
    cdfFilename = [metadata.filename '.cdf'];
end

disp(['Writing ancillary sensor data to the netCDF file ' cdfFilename])

defineBurstVECcdfFile(cdfFilename,metadata);

%write the sensor data first
ncwrite(cdfFilename,'lat',metadata.latitude);
ncwrite(cdfFilename,'lon',metadata.longitude);
ncwrite(cdfFilename,'time',time);
ncwrite(cdfFilename,'time2',time2);
ncwrite(cdfFilename,'Temperature',SEN.temperature);
ncwrite(cdfFilename,'Battery',SEN.battery);
ncwrite(cdfFilename,'Pitch',SEN.pitch);
ncwrite(cdfFilename,'Roll',SEN.roll);
ncwrite(cdfFilename,'Heading',SEN.heading);


%now get the burst data from the .dat file and write incrementally   
n = instmeta.VECSamplesPerBurst;
datid = fopen([baseFile '.dat']);
str = fgetl(datid);buf = textscan(str,'%f');bstnum = buf{1}(1);
%determine number of bytes per line of data for rewinding purposes
byts = ftell(datid);
%advance to the first good burst
while bstnum<goodbursts(1)
    str = fgetl(datid);buf = textscan(str,'%f');bstnum = buf{1}(1);
end
%now go back one line, because the above carries us one line too many
fseek(datid,byts*-1,'cof');
%we can't assume that each burst has the right number of datapoints, so
%we'll set up a structure that does have the right number, and fill them in
%line by line, and any missing samples will get a fill values
% fields in the output structure
RAW = struct('BURST',NaN(n,1),'ENS',NaN(n,1),'V1',NaN(n,1),'V2',NaN(n,1),...
    'V3',NaN(n,1),'AMP1',NaN(n,1),'AMP2',NaN(n,1),'AMP3',NaN(n,1),...
    'SNR1',NaN(n,1),'SNR2',NaN(n,1),'SNR3',NaN(n,1),'COR1',NaN(n,1),...
    'COR2',NaN(n,1),'COR3',NaN(n,1),'PRESS',NaN(n,1),'ANALOG1',NaN(n,1),...
    'ANALOG2',NaN(n,1),'CHK',NaN(n,1));
flds = fieldnames(RAW);
tic
disp('Now reading in burst data')
for i = 1:nbursts
    count = 1;
    while bstnum == goodbursts(i);
        str = fgetl(datid);buf = textscan(str,'%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f');
        ID = buf{2};
        for j = 1:length(flds)
            fld = flds{j};
%             RAW.(fld)(count,i) = buf{j};
            RAW.(fld)(ID,i) = buf{j};
        end
        bstnum = RAW.BURST(count,i);
        if count == spb
            RAW.Burst(i) = bstnum;
            bstnum = bstnum+1;
            RAW.SPB(i) = count;
            continue
        else
            count = count +1;
        end
    end
        
% TO simplify, I'll save the writing of the data until after all data have 
%   been read in. This allows the user to save the raw data as a structure
%   as well.
%     %convert from m/s to cm/s
%     A.VEL1 = A.V1.*100;
%     A.VEL2 = A.V2.*100;
%     A.VEL3 = A.V3.*100;
%     
%     ncwrite(cdfFilename,'VEL1',A.VEL1(:,i),[1 i]);
%     ncwrite(cdfFilename,'VEL2',A.VEL2(:,i),[1 i]);
%     ncwrite(cdfFilename,'VEL3',A.VEL3(:,i),[1 i]);
%     ncwrite(cdfFilename,'AMP1',A.AMP1(:,i),[1 i]);
%     ncwrite(cdfFilename,'AMP2',A.AMP2(:,i),[1 i]);
%     ncwrite(cdfFilename,'AMP3',A.AMP3(:,i),[1 i]);
%     ncwrite(cdfFilename,'Pressure',A.PRESS(:,i),[1 i]);
%     if isfield(metadata,'AnalogInput1')
%         ncwrite(cdfFilename,'AnalogInput1',A.ANALOG1(:,i),[1 i]);
%     end
%     if isfield(metadata,'AnalogInput2')
%         ncwrite(cdfFilename,'AnalogInput2',A.ANALOG2(:,i),[1 i]);
%     end

    if rem(i,100)==0
        disp(['Finished reading in ',num2str(i),' bursts'])
    end
end

%convert to cm/s according to scaling in the status code - not sure this is
%correct so I will leave it out for now, and assume all data start in m/s
RAW.VEL1 = RAW.V1.*100;RAW.VEL2 = RAW.V2.*100;RAW.VEL3 = RAW.V3.*100;
% switch scaling
%     case '0'
%         RAW.VEL1 = RAW.V1.*100;RAW.VEL2 = RAW.V2.*100;RAW.VEL3 = RAW.V3.*100;
%     case '1'
%         RAW.VEL1 = RAW.V1.*1000;RAW.VEL2 = RAW.V2.*1000;RAW.VEL3 = RAW.V3.*1000;
% end
RAW.time = time;RAW.time2 = time2;
RAW.jd = RAW.time + RAW.time2/86400000;

% now do the fill values for the RAW data
N = length(RAW.time);flds = fieldnames(RAW);
for i = 1:length(RAW)
    fld = flds{i};
    if length(RAW.(fld) == N); %i.e. only work on timeseries
        if strcmp(fld,{'time';'time2'}),
            continue
        else
            nanind = isnan(RAW.(fld));RAW.(fld)(nanind) = doublefill;
        end
    end
end

% tmp(1,1,:,:) = RAW.VEL1;ncwrite(cdfFilename,'VEL1',tmp);
% tmp(1,1,:,:) = RAW.VEL2;ncwrite(cdfFilename,'VEL2',tmp);
% tmp(1,1,:,:) = RAW.VEL3;ncwrite(cdfFilename,'VEL3',tmp);
% tmp(1,1,:,:) = RAW.AMP1;ncwrite(cdfFilename,'AMP1',tmp);
% tmp(1,1,:,:) = RAW.AMP2;ncwrite(cdfFilename,'AMP2',tmp);
% tmp(1,1,:,:) = RAW.AMP3;ncwrite(cdfFilename,'AMP3',tmp);
% tmp(1,1,:,:) = RAW.COR1;ncwrite(cdfFilename,'COR1',tmp);
% tmp(1,1,:,:) = RAW.COR2;ncwrite(cdfFilename,'COR2',tmp);
% tmp(1,1,:,:) = RAW.COR3;ncwrite(cdfFilename,'COR3',tmp);
% tmp(1,1,:,:) = RAW.PRESS;ncwrite(cdfFilename,'Pressure',tmp);

% get the sensor depth and height figured out so there is no confusion
% Per email with EM, 7/17/2013, it was decided that any instrument that
% doesn't have all sensors co-located will not get a global attribute of
% 'initial_instrument_height' only the sensors will get heights and depths,
% but there is still the coordinate variable of 'Depth' so it was decided
% that would be the depth of the velocity measurement, which includes the
% semple volume offset;
% first check to make sure we have the height of the pressure sensor -
% If not, assume a fixed probe, and determine the height;
if ~isfield(metadata,'probe_height')
    disp('No probe height found in the metadata')
    metadata.probe_height = input('Enter height of the Vector probe, center transducer: ');
    
end
    
if ~isfield(metadata,'pressure_sensor_height')
    switch metadata.orientation
        case 'UP'           %the probe is facing up, pressure sensor below it
            metadata.pressure_sensor_height = metadata.probe_height - 0.215;
            disp('No pressure sensor height in the metadata, thus it is assumed to be a fixed vector')
            disp('---The orientation is UP, thus the pressure sensor height is set to the probe height - 0.215m')
        case 'DOWN'        %the probe is facing up, pressure sensor below it
            metadata.pressure_sensor_height = metadata.probe_height + 0.215;
            disp('No pressure sensor height in the metadata, thus it is assumed to be a fixed vector')
            disp('---The orientation is DOWN, thus the pressure sensor height is set to the probe height + 0.215m')
    end
end
% now get the height of the velocity measurement
% need to double check orientation; CMG convention of orientation is that
% when the head is pointing up, the instrument is 'UP', and vice versa. The
% orientation is stored in the status code upon each 'new' initialization,
% which for bursts means each burst, but in continuous mode, it is only
% when it first starts up. We also need to swap the sign of the Y and Z
% directions if instrument was in XYZ coordinates, and the orientation of 
% the probe is opposite to the orientation of the compass module.

% we only need to flip the transformation matrix if the instrument is
% pointing down, and the matrix is only used to go from beam coordinates to
% XYZ (or ENU).

T = instmeta.VECTransformationMatrix;


switch metadata.orientation
    case 'UP'           %the probe is facing up, pressure sensor below it
        switch orient
            case '0'
                disp('User instructed that instrument was pointing UP, instrument states it was DOWN; using user input for tranforms')
            case '1'
                disp('User instructed that instrument was pointing UP, instrument states it was UP')
        end
        metadata.velocity_height = metadata.probe_height + 0.157;
    case 'DOWN'        %the probe is facing down, pressure sensor above it
        switch orient
            case '0'
                disp('User instructed that instrument was pointing DOWN, instrument states it was DOWN')
                disp('Flipping Transformation Matrix')
                if strcmp(instmeta.VECCoordinateSystem,'XYZ')
                    disp('Coordinate system is XYZ and user states that probe is down - changing sign of Y and Z directions')
                    RAW.VEL2 = RAW.VEL2*-1;RAW.VEL3 = RAW.VEL3*-1;
                end
                T(2,:) = -T(2,:);
                T(3,:) = -T(3,:);
            case '1'
                disp('User instructed that instrument was pointing DOWN, instrument states it was UP; using user input for tranforms')
                disp('Flipping Transformation Matrix')
                T(2,:) = -T(2,:);
                T(3,:) = -T(3,:);
                if strcmp(instmeta.VECCoordinateSystem,'XYZ')
                    RAW.VEL2 = RAW.VEL2*-1;RAW.VEL3 = RAW.VEL3*-1;
                    disp('Coordinate system is XYZ and user states that probe is down - changing sign of Y and Z directions')
                end
        end
        metadata.velocity_height = metadata.probe_height - 0.157;
end
ncwrite(cdfFilename,'burst',RAW.Burst);
ncwrite(cdfFilename,'VEL1',RAW.VEL1);
ncwrite(cdfFilename,'VEL2',RAW.VEL2);
ncwrite(cdfFilename,'VEL3',RAW.VEL3);
ncwrite(cdfFilename,'AMP1',RAW.AMP1);
ncwrite(cdfFilename,'AMP2',RAW.AMP2);
ncwrite(cdfFilename,'AMP3',RAW.AMP3);
ncwrite(cdfFilename,'COR1',RAW.COR1);
ncwrite(cdfFilename,'COR2',RAW.COR2);
ncwrite(cdfFilename,'COR3',RAW.COR3);
ncwrite(cdfFilename,'Pressure',RAW.PRESS);
ncwrite(cdfFilename,'TransMatrix',T);
% convert the analog input data to volts before writing
if isfield(metadata,'AnalogInput1')
    RAW.analog1 = RAW.ANALOG1*5/65536;
    ncwrite(cdfFilename,'AnalogInput1',RAW.analog1);
end
if isfield(metadata,'AnalogInput2')
    RAW.analog2 = RAW.ANALOG2*5/65536;
    ncwrite(cdfFilename,'AnalogInput2',RAW.analog1);
end

disp(['Finished reading and writing ',num2str(nbursts),' bursts in ',num2str(toc/60),' minutes'])


Depth = nanmean(nanmean(RAW.PRESS));
wdepth = Depth + metadata.pressure_sensor_height;
metadata.WATER_DEPTH = wdepth;
metadata.WATER_DEPTH_source = 'water depth = MSL from pressure sensor';
metadata.WATER_DEPTH_datum = 'MSL';
RAW.Depth = metadata.WATER_DEPTH - metadata.velocity_height;
        
%now fix some of the attributes based on the data just read in
ncwriteatt(cdfFilename,'/','WATER_DEPTH',metadata.WATER_DEPTH);
ncwriteatt(cdfFilename,'/','WATER_DEPTH_source',metadata.WATER_DEPTH_source);
ncwriteatt(cdfFilename,'/','WATER_DEPTH_datum',metadata.WATER_DEPTH_datum);
ncwrite(cdfFilename,'depth',RAW.Depth);
ncwriteatt(cdfFilename,'depth','initial_instrument_height',metadata.velocity_height);
ncwriteatt(cdfFilename,'depth','nominal_instrument_depth',RAW.Depth);


%% fix some of the attributes
% fix the sensor depth and height attributes - don't want to apply to every
% variable, as the analog sensors have their own, and the probe may be a 
% different height than the electronics and hence pressure and temperature
% so need a list of variables to apply each to
velsens = {'VEL1';'VEL2';'VEL3';'AMP1';'AMP2';'AMP3';'COR1';'COR2';'COR3'};
for i = 1:length(velsens)
    fld = velsens{i};
    ncwriteatt(cdfFilename,fld,'initial_sensor_height',metadata.velocity_height);
    ncwriteatt(cdfFilename,fld,'nominal_sensor_depth',RAW.Depth);
end
ancsens = {'Pitch';'Roll';'Heading';'Temperature';'Pressure'};
for i = 1:length(ancsens)
    fld = velsens{i};
    ncwriteatt(cdfFilename,fld,'initial_sensor_height',metadata.pressure_sensor_height);
    ncwriteatt(cdfFilename,fld,'nominal_sensor_depth',metadata.WATER_DEPTH-metadata.pressure_sensor_height);
end

%% define the file


%% --------------------------------