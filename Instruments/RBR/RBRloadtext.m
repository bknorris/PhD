%RBR Solo/Duet Load n Save
clear %a good idea

%since the Duet and Pressure loggers couldn't be exported as .mat, I have
%the fun opportunity to read in the data file and compile all the metadata.
%This code will read the textfile line by line and get the relevant
%information and data.

%define metadata
inst_type = 'Solo83P'; %inst type for saving data
metadata.serial = 77683;
metadata.inst_type = 'RBR Solo Pressure Gauge';  % type of instrument and instrument manufacturer
metadata.data_cmt = 'Mekong 2014 Dense Pneumatophore Study Solo 83P'; % any comment
metadata.latitude = 9.565555556;
metadata.longitude = 106.2923889;
metadata.hab = 0; %instrument height above bottom (in mm)
metadata.deployment_date = [2014 10 02 18 00 00]; %times to crop deployment
metadata.recovery_date = [2014 10 03 12 00 00];
file = '077683_20141003_1633.txt'; %text file to be loaded

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fname = fopen(file);
str = fgetl(fname);
hdr = cell(30,1); %estimate number of lines of header to be 30
count = 1;
%read the header
while ~strncmp(str,'NumberOfSamples',15)
    str = fgetl(fname);
    if isempty(str)
        str = fgetl(fname);
    end
    hdr(count,1) = textscan(str,'%s','delimiter','\n');
    if strncmp(str,'LoggingStartTime',16)
        is=strfind(str,'=');
        metadata.starttime = str(is+1:end);
    elseif strncmp(str,'LoggingEndTime',14)
        is=strfind(str,'=');
        metadata.endtime = str(is+1:end);
    elseif strncmp(str,'LoggingSamplingPeriod',21)
        is=strfind(str,'=');
        metadata.samprate = str(is+1:end);
    elseif strncmp(str,'NumberOfChannels',16)
        is=strfind(str,'=');
        metadata.num_chan = str2double(str(is+1:end));
    end
    count = count+1;
end
hdr = [hdr{:}];
str = fgetl(fname); %skip to next line
hdr = hdr(~cellfun('isempty',hdr)); %remove empty cells

%get channel and event metadata
chan = strncmp(hdr,'Channel',7)';chanid = find(chan == 1);
for i = 1:length(chanid)
    Chan = char(hdr(chanid(i)));
    flds = regexp(Chan,'\.(\w*)\=','tokens');
    flds = char(flds{:});
    channum = regexp(Chan,'\[(\d*)\]','tokens');
    eqls = strfind(Chan,'=');
    chanval = Chan(eqls+1:end);
    strct = char(strcat('channel',channum{:}));
    channels.(strct).(flds) = chanval;
end
event = strncmp(hdr,'Event',5)';eventid = find(event == 1);
for i = 1:length(eventid)
    Event = char(hdr(eventid(i)));
    flds = regexp(Event,'\.(\w*)\=','tokens');
    flds = char(flds{:});
    eventnum = regexp(Event,'\[(\d*)\]','tokens');
    eqls = strfind(Event,'=');
    evntval = Event(eqls+1:end);
    strct = char(strcat('event',eventnum{:}));
    events.(strct).(flds) = evntval;
end
clearvars Chan Event flds eqls chanval eventval strct
metadata.channels = channels;
metadata.events = events;
sampn = strfind(hdr{end},'=');
metadata.nsamp = str2double(hdr{end}(sampn+1:end));

%now, get the fieldnames of the channels
str = fgetl(fname);chanflds = regexp(str, '  +', 'split');
chanflds(strcmp('',chanflds)) = [];chanflds = strrep(chanflds,' ',''); %remove any whitespace
chanflds = regexprep(chanflds,'&',''); %fix date & time string
%replace given variable names
chanflds = strrep(chanflds,'DateTime','Datetime');
chanflds = strrep(chanflds,'temp08','Temp');
chanflds = strrep(chanflds,'pres03','Pres');
chanflds = strrep(chanflds,'pres08','SeaPres');
chanflds = strrep(chanflds,'pres16','Pres');
chanflds = strrep(chanflds,'dpth01','Depth');
%define start and end times
deploy = datenum(metadata.deployment_date);deploy = datestr(deploy,'dd-mmm-yyyy HH:MM:SS.FFF');
recover = datenum(metadata.recovery_date);recover = datestr(recover,'dd-mmm-yyyy HH:MM:SS.FFF');

%now read in the data
disp('Now reading the data...')
tic
n = metadata.nsamp;
datafill = NaN(n,1);
for i = 1:length(chanflds)
    DAT.(chanflds{i}) = datafill; %preallocate structure DAT
end
flds = fieldnames(DAT);
FRMT = repmat('%f32',1,length(flds)-1);FRMT = ['%*s%*s' FRMT];
dateFRMT = repmat('%*n',1,length(flds)-1);dateFRMT = ['%s%s' dateFRMT]; %use this to extract the dates only

%move file counter to the line that matches the specified deployment time
startid = 0;
while startid < 1;
    str = fgetl(fname);
    if strcmp(str(1:24),deploy)
        startid = 1;
    end
end
%read in the data line-by-line until reaching the line where the date
%matches the specified recovery time
endid = 0;
i = 1;
while endid < 1;
    %in case there are null values, add a catch and fill with Inf
    if strfind(str,'null')
        nullid = strfind(str,'null');
        str(nullid:nullid+3) = '1E99';
    end
    %since the formatting for these text files is a bit strange, we will
    %have to use a workaround: two textscan calls, one with formatting only
    %for the datetime, and the other with formatting only for the rest of
    %the data.
    dates = textscan(str,dateFRMT,'whitespace','  ');
    datetime = strcat(dates{1},{' '},dates{2}); %concatenate with whitespace
    dvec = datevec(datetime,'dd-mmm-yyyy HH:MM:SS.FFF'); %vector for yearday
    datetime = datenum(datetime,'dd-mmm-yyyy HH:MM:SS.FFF'); %serial datenum for datetime
    dat = textscan(str,FRMT,'whitespace','  '); %read the data
    buf = [datetime dat];
    yy = dvec(:,1);mo = dvec(:,2);dd = dvec(:,3); %create yearday vector for Julia
    hh = dvec(:,4);mi = dvec(:,5);ss = dvec(:,6);
    DAT.Yearday(i,:) = yearday(yy,mo,dd+hh/24+mi/(24*60)+ss/(24*60*60));
    for j = 1:length(flds)
        fld = flds{j};
        
        DAT.(fld)(i,:) = buf{j};
    end
    if strcmp(str(1:24),recover)
        endid = 1;
        disp(['Finished reading in ',num2str(i),' samples in ',num2str(toc/60),' minutes'])
    end
    str = fgetl(fname); %get next line
    if (i/10000)-floor(i/10000)==0 %count samples read
        disp([num2str(i) ' samples read']),
    end
    i = i+1; %count up to next index
end
%remove any extra leading and/or trailing NaNs, create final structure 
for i = 1:length(flds)
    nans = isnan(DAT.(flds{i}));
    DAT.(flds{i})(nans) = [];
end
RBR = DAT;
RBR.Metadata = metadata;
clearvars -except RBR recover inst_type
sname = [inst_type '_' datestr(recover,'ddmmyy')];
save(sname,'RBR','-v7.3')
disp(['File saved as ' sname '.mat'])

