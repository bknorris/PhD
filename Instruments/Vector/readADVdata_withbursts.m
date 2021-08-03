function adv = readADVdata(filename,metadata)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%A function to read in raw Nortek Vector single point current meter data
%and create a structure (adv) which contains the metadata, .sen file data,
%and the .dat file data.
%
%   filename - the base filename without the terminal string
%   metadata - a structure of metadata, user defined. See adv_dataprocess.m
%   for an example of the required inputs
%
%  Contains an adaptation of vec2cdfBurstNative.m
%  Written by Kurt J Rosenberger
%  krosenberger@usgs.gov
%  USGS Pacific Coastal Marine Science Center
%
%  This script was compiled by Benjamin K Norris, 2015
%  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if nargin<1
    help(mfilename)
    return
end
thefillvalue = NaN; %define a fill value for non-elements
instmeta = readVEChdr(filename);

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
sendata = load(senFile);[n,m] = size(sendata);
%Find breaks in the timebase to delineate bursts. The sensor data is
%collected at 1Hz, regardless of the velocity sampling rate.
yy = instmeta.senflds.Year;mo = instmeta.senflds.Month;
dd = instmeta.senflds.Day;hh = instmeta.senflds.Hour;
mi = instmeta.senflds.Minute;ss = instmeta.senflds.Second;
dn = datenum([sendata(:,yy) sendata(:,mo) sendata(:,dd) sendata(:,hh) sendata(:,mi) sendata(:,ss)]);
ddn = find(diff(dn)>(10/86400))+1; %finds the start indices of all the bursts
bstdn = dn(ddn);
gooddata = find(dn>=deploy&dn<=recover); %indices of good data

if strcmp(instmeta.burstint,'CONTINUOUS') 
    %if mode is continuous, there are no breaks for bursts. This means that
    %gooddata should be used to crop datafiles to start and stop times as
    %defined as metadata.deployment_date and metadata.recovery_date.
    startind = gooddata;
    nsamp = length(gooddata);
    datadn = dn(gooddata);
    spb = instmeta.spb;
    disp(['There are ',num2str(nsamp),' good samples in the file'])
else
    goodbursts = find(bstdn>=deploy&bstdn<=recover);
    burstdn = bstdn(goodbursts);
    startind = ddn(goodbursts);
    %check to see if the instrument was stopped mid-burst
    spb = instmeta.spb;
    if (gooddata(end)-startind(end))<spb
        goodbursts = goodbursts(1:end-1);burstdn = burstdn(1:end-1);startind = startind(1:end-1);
    end
    if isfield(metadata,'good_bursts')
        disp(['User instructed to process bursts ',num2str(metadata.good_bursts(1)),' to ',num2str(metadata.good_bursts(2))])
        ind1 = near(goodbursts,metadata.good_bursts(1),1);
        ind2 = near(goodbursts,metadata.good_bursts(2),1);
        burstdn = burstdn(ind1:ind2);
        startind = startind(ind1:ind2);
        goodbursts = metadata.good_bursts(1):1:metadata.good_bursts(2);
        nsamp = length(goodbursts);
    else
        nsamp = length(goodbursts);
        disp(['There are ',num2str(nsamp),' good bursts in the file'])
    end
    gooddata = goodbursts; %save same variable names to avoid conflict later
    datadn = burstdn;
end

%get the error and status code from the first good measurement to determine
%orientation and status information for the instrument
senfid = fopen(senFile);
for i = 1:startind(1)+1
    str = fgetl(senfid);
end
X = regexp(str,' +','split');
error = char(X{7});status = char(X{8});
orient = (status(8));scaling = (status(7));
fclose(senfid);

%Crop sensor data to gooddata, load into a structure
senfields = fieldnames(instmeta.senflds);
SENdata = sendata(1:nsamp,:);
for i = 8:16
    senfld = senfields{i};
    SEN.(senfld) = SENdata(:,i);
end
SEN.datenum = datadn;
N = length(SEN.datenum);fields = fieldnames(SEN);
SEN.int = instmeta.burstinterval;
% put in fill values    
for i = 1:length(SEN)
    fld = fields{i};
    if length(SEN.(fld) == N); %i.e. only work on timeseries
        if strcmp(fld,{'time';'time2'}),
            continue
        else
            nanind = isnan(SEN.(fld));SEN.(fld)(nanind) = thefillvalue;
        end
    end
end

SEN.time = floor(datadn); SEN.time2 = (datadn - floor(datadn))*86400000;
%define metadata fields
metadata.sampleinds = gooddata;
metadata.inst_type = 'Nortek Vector Acoustic Doppler Velocimeter';
metadata.instmeta = instmeta;
metadata.nsamp = nsamp;
metadata.samplerate = SEN.int;
metadata.burstlength = spb;
metadata.start_time = datestr(SEN.datenum(1),'yyyy-mm-dd HH:MM:SS:FFF');
metadata.stop_time = datestr(SEN.datenum(end),'yyyy-mm-dd HH:MM:SS:FFF');

%create new burst timebase
time = zeros(spb,nbursts);time2 = zeros(spb,nbursts);
dTb = 1/instmeta.samprate*1000;
bsttime = 0:dTb:dTb*(spb-1);
for i = 1:nbursts
    time(:,i) = repmat(SEN.time(i),1,spb);
    time2(:,i) = SEN.time2(i) + bsttime;
end
% sen = fopen(senFile);
% count = 1;
% while ~feof(sen)
%     str = fgetl(sen);
%     X = regexp(str,' +','split');
%     day(count,:) = str2num(X{SENFIELDS.Day});
%     month(count,:) = str2num(X{SENFIELDS.Month});
%     year(count,:) = str2num(X{SENFIELDS.Year});
%     hour(count,:) = str2num(X{SENFIELDS.Hour});
%     minute(count,:) = str2num(X{SENFIELDS.Minute});
%     second(count,:) = str2num(X{SENFIELDS.Second});
%     metadata.heading(count,:) = str2num(X{SENFIELDS.Heading});    
%     metadata.pitch(count,:) = str2num(X{SENFIELDS.Pitch});
%     metadata.roll(count,:) = str2num(X{SENFIELDS.Roll});
%     metadata.batt_voltage(count,:) = str2num(X{SENFIELDS.Battvolt});
%     count = count+1;
% end

%write timestamps from .sen file to the adv structure
% adv.datenum = datenum(year,month,day,hour,minute,second);
% adv.yearday = yearday(year,month,day+hour/24+minute/(24*60)+second/(24*60*60)); %to make Julia happy
%These are 1Hz, we need timestamps to match the sample rate
end
