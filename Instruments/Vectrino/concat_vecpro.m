function VPRO = concat_vecpro(fileList,newFile)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%concat_vecpro concatenates several Nortek Vectrino II Profiler files (in
%.mat) and deals appropriately with time gaps. % files must have the same
%variables, attributes, i.e. be a part of the same instrument time series.
%
% Usage: concatcdf_gaps(fileList, newFile)
% Where:
%       fileList is a cell array of file names (use dir)
%       newFile is the name of the resulting file
%
% Calls cmgbridge.m,cmgidgaps.m,cmgdataclean.m
%
% Contains a partial adaptation of concatcdf_gaps.m
% Written by Kurt J Rosenberger
% krosenberger@usgs.gov
% USGS Pacific Coastal Marine Science Center
%
% This script was written by Benjamin K Norris, 2015
% University of Waikato, New Zealand
%
% UPDATES:
% 08/07/2015: added bottom check timebase and bottom check to fields
% bridged by the program.
% 07/10/2015: fixed issue with concatenating data files; timebase and data
% were different dimensions due to an error in grouping data files. All
% other concatenation routines have been updated to this method.
% 16/10/2015: added 'num' variable to the first for loop that saves the
% length of the files to be loaded. 'id' is the index of the longest variable
% and is used to preallocate the structures into which data is loaded.
% 19/10/2015: added a check for the bottom distance routine to 'try'
% loading data from the first file. There are some VecPro files where the
% bottom distance subfields are not written to the matlab structure, so
% this routine checks if these fields exist and save out the BD time series
% as well as the BD field.
% 03/02/2016: For instruments with long gaps in between sampling periods,
% I've added a catch that looks for the last time of the current file and
% compares it with the first time of the next file, then fills the gaps
% with NaNs accordingly. NaNs are bridged linearly.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[path,name,ext] = fileparts(fileList{1});

if nargin < 2
    disp('concat_vecpro requires two inputs')
    help(mfilename);
    return
end
if ~exist('fileList', 'var')
    disp('A list of two or more files is required')
    return
end
if ~exist('newFile', 'var')
    newFile = inputdlg('Enter a file name');fname = newFile{:};
else
    fname = newFile;
end

if length(fileList)<2, disp('Need more than one file...'); return; end

%get the data
disp('Reading in the data from multiple files')
for ifile = 1:length(fileList),
    A(ifile) = load(fileList{ifile});
    num(ifile) = length(A(ifile).Data.Profiles_HostTimeMatlab);
end
%get the maximum individual file size to be read
[~,id] = max(num);

%open one of the original files to get variable names
a = load(fileList{1});
datanames = fieldnames(a.Data);
clearvars a
%first get the times from all the files to create a new timebase
time = A(id).Data.Profiles_HostTimeMatlab'; %time variable needs to be as large as the largest file being read
timev = NaN(length(time)*length(A),1);
%UPDATE: 03/02/2016 I need to add in a catch for timebase issues. Some
%VecPros have timebases that jump due to powerloss, etc. t1 and t2 will
%designate the first and last times of the subsequent and current file
%being concatenated, respectively.
fs = A(1).Config.sampleRate; %in Hz (samples/second)
tgap = datenum(0,0,0,0,0,2); %allowable gap in time: 2 seconds of data (for 50Hz instruments that's 100 samples)
for i = 1:length(A)
    [~,n] = size(A(i).Data.Profiles_HostTimeMatlab);
    if i == 1
        idx = 1:n;
    elseif i > 1
        t1 = A(i).Data.Profiles_HostTimeMatlab(1);
        %check gaps in timebase
        if diff([t2 t1]) > tgap
            %if tgap is larger than 2 sec, calc num of samples in the gap
            nsec = etime(datevec(t1),datevec(t2));
            %b/c of rounding, some error in timebase will occur. This is on the order of < 1/fs, so I'll consider it acceptable. 
            nsamp = round(nsec*fs);
            idx = nanidx+nsamp:nanidx+n+nsamp; %add the number of gaps as NaNs to the first row after the last file
        elseif diff([t2 t1]) < tgap
            idx = nanidx:nanidx+n; %new index is NaNindex to the length of the next data file
        end
    end
    for ii = 1:length(A(i).Data.Profiles_HostTimeMatlab');
        timev(idx(ii),:) = A(i).Data.Profiles_HostTimeMatlab(:,ii);
    end
    nanidx = idx(end);
    t2 = A(i).Data.Profiles_HostTimeMatlab(end);
end

%bridge gaps in the time vector with cmgbridge
nlin = 1E6;
maxgaps = 1E8; %needs to be large to catch files with many gaps
[time,nogap] = cmgbridge(timev,maxgaps,1,maxgaps); %fill gaps linearly so the next iteration identifies the next set of NaNs correctly

time(isnan(time(:,1)),:) = [];
if nogap
    disp('Gaps in time record filled')
end
clearvars newdata data n m idx

%now loop through all the record variables and paste the data and fill gaps
disp('Concatenating all the data and filling gaps with cmgbridge')
DATA = struct();
DATA.Time = time; %write timebase to file
fields = fieldnames(A(1).Data);
nums = 1:27;
for j = nums %only fill gaps in the data fields
    ind = find(strcmp(datanames(j),fields));
    [p,q] = size(A(id).Data.(fields{ind}));
    if p == 1
        data = NaN(length(time),p);
        for i = 1:length(A)
            [~,n] = size(A(i).Data.(fields{ind}));
            if i == 1
                idx = 1:n;
            elseif i > 1
                %need to find the timegaps again for the data concatenation
                %step
                t1 = A(i).Data.Profiles_HostTimeMatlab(1);
                if diff([t2 t1]) > tgap
                    %if tgap is larger than 2 sec, calc num of samples in the gap
                    nsec = etime(datevec(t1),datevec(t2));
                    %b/c of rounding, some error in timebase will occur. This is on the order of < 1/fs, so I'll consider it acceptable.
                    nsamp = round(nsec*fs);
                    idx = nanidx+nsamp:nanidx+n+nsamp; %add the number of gaps as NaNs to the first row after the last file
                elseif diff([t2 t1]) < tgap
                    idx = nanidx:nanidx+n; %new index is NaNindex to the length of the next data file
                end
            end
            for ii = 1:length(A(i).Data.(fields{ind}))
                data(idx(ii),:) = A(i).Data.(fields{ind})(:,ii);
            end
            nanidx = idx(end);
            t2 = A(i).Data.Profiles_HostTimeMatlab(end);
        end
        %bridge gaps with cmgbridge
        [newdata,nogap] = cmgbridge(data,nlin,nlin,maxgaps);
        newdata(isnan(newdata(:,1)),:) = [];
        data = newdata;
    elseif p > 1
        data = NaN(length(time),q);
        for i = 1:length(A)
            [n,~] = size(A(i).Data.(fields{ind}));
            if i == 1
                idx = 1:n;
            elseif i > 1
                t1 = A(i).Data.Profiles_HostTimeMatlab(1);
                if diff([t2 t1]) > tgap
                    %if tgap is larger than 2 sec, calc num of samples in the gap
                    nsec = etime(datevec(t1),datevec(t2));
                    %b/c of rounding, some error in timebase will occur. This is on the order of < 1/fs, so I'll consider it acceptable.
                    nsamp = round(nsec*fs);
                    idx = nanidx+nsamp:nanidx+n+nsamp; %add the number of gaps as NaNs to the first row after the last file
                elseif diff([t2 t1]) < tgap
                    idx = nanidx:nanidx+n; %new index is NaNindex to the length of the next data file
                end
            end
            for ii = 1:length(A(i).Data.(fields{ind}))
                data(idx(ii),:) = A(i).Data.(fields{ind})(ii,:);
            end
            nanidx = idx(end);
            t2 = A(i).Data.Profiles_HostTimeMatlab(end);
        end
        %leave NaNs as gaps so future data analysis doesn't include gaps as
        %values
        data(nanidx:length(data),:) = []; %clips trailing NaNs
    end
    %identify gaps for user
    if cmgidgaps(data(:,1)) > 0
        [gaps,~,~,lgap] = cmgidgaps(data(:,1));
        disp(['There are ' num2str(gaps) ' gaps in ' fields{ind}])
        disp(['Gaps are ' num2str(lgap) ' samples long'])
    end
    DATA.(fields{ind}) = data;
    clearvars newdata data n m idx
end

%For scripts that require the bottom distance, concatenate these fields
%together using the same routine (this is not done in the step above).
%need to get the bottom check timebase and bottom check

%UPDATE 19/10/2015: some files being loaded do not have the bottom distance
%subfields in the Data structure. I've inserted a 'try' command to see if
%this field exists in the largest input array. If it doesn't exist, then the
%script searches for the next largest input array until a bottom distance
%subfield is found.
bcti = find(strcmp(fields,'BottomCheck_HostTimeMatlab'));
bci = find(strcmp(fields,'BottomCheck_BottomDistance'));
try
    bct = A(id).Data.BottomCheck_HostTimeMatlab';
    bctv = NaN(length(bct)*2*length(A),1);
catch
    count = 1;
    bct = [];
    while isempty(bct)
        if count == 1
            num(id) = 0;
            [~,newid] = max(num);
            if isfield(A(newid).Data,'BottomCheck_HostTimeMatlab')
                bct = A(newid).Data.BottomCheck_HostTimeMatlab';
                bctv = NaN(length(bct)*length(A),1);
            end
        elseif count > 1
            num(newid) = 0;
            [~,newid] = max(num);
            if isfield(A(newid).Data,'BottomCheck_HostTimeMatlab')
                bct = A(newid).Data.BottomCheck_HostTimeMatlab';
                bctv = NaN(length(bct)*length(A),1);
            end
        end
        count = count+1;
    end
end
for i = 1:length(A)
    if isfield(A(i).Data,'BottomCheck_HostTimeMatlab')
        [~,n] = size(A(i).Data.BottomCheck_HostTimeMatlab);
        if i == 1
            idx = 1:n;
        elseif i > 1
            t1 = A(i).Data.Profiles_HostTimeMatlab(1);
            if diff([t2 t1]) > tgap
                %if tgap is larger than 2 sec, calc num of samples in the gap
                nsec = etime(datevec(t1),datevec(t2));
                %note the bottom check only writes at 10 Hz
                nsamp = round(nsec*(10));
                idx = nanidx+nsamp:nanidx+n+nsamp; %add the number of gaps as NaNs to the first row after the last file
            elseif diff([t2 t1]) < tgap
                idx = nanidx:nanidx+n; %new index is NaNindex to the length of the next data file
            end
        end
        for ii = 1:length(A(i).Data.BottomCheck_HostTimeMatlab');
            bctv(idx(ii),:) = A(i).Data.BottomCheck_HostTimeMatlab(:,ii);
        end
        nanidx = idx(end);
        t2 = A(i).Data.Profiles_HostTimeMatlab(end);
    else
        continue
    end
end
%bridge gaps with cmgbridge
[newdata,nogap] = cmgbridge(bctv,nlin,maxgaps,maxgaps);
newdata(isnan(newdata(:,1)),:) = [];
if ~nogap
    disp(['Gaps in ' fields{bcti} ' filled'])
end
DATA.BottomCheck_HostTimeMatlab = newdata;
clearvars newdata data n m idx

try
    bc = A(id).Data.BottomCheck_BottomDistance;
    bcv = NaN(length(bc)*2*length(A),1);
catch
    count = 1;
    bc = [];
    while isempty(bc)
        if count == 1
            num(id) = 0;
            [~,newid] = max(num);
            if isfield(A(newid).Data,'BottomCheck_BottomDistance')
                bc = A(newid).Data.BottomCheck_BottomDistance';
                bcv = NaN(length(bc)*length(A),1);
            end
        elseif count > 1
            num(newid) = 0;
            [~,newid] = max(num);
            if isfield(A(newid).Data,'BottomCheck_BottomDistance')
                bc = A(newid).Data.BottomCheck_BottomDistance';
                bcv = NaN(length(bc)*length(A),1);
            end
        end
        count = count+1;
    end
end

for i = 1:length(A)
    if isfield(A(i).Data,'BottomCheck_BottomDistance');
        [n,~] = size(A(i).Data.BottomCheck_BottomDistance);
        if i == 1
            idx = 1:n;
        elseif i > 1
            t1 = A(i).Data.Profiles_HostTimeMatlab(1);
            if diff([t2 t1]) > tgap
                %if tgap is larger than 2 sec, calc num of samples in the gap
                nsec = etime(datevec(t1),datevec(t2));
                nsamp = round(nsec*10);
                idx = nanidx+nsamp:nanidx+n+nsamp; %add the number of gaps as NaNs to the first row after the last file
            elseif diff([t2 t1]) < tgap
                idx = nanidx:nanidx+n; %new index is NaNindex to the length of the next data file
            end
        end
        for ii = 1:length(A(i).Data.BottomCheck_BottomDistance');
            bcv(idx(ii),:) = A(i).Data.BottomCheck_BottomDistance(ii,:);
        end
        nanidx = idx(end);
        t2 = A(i).Data.Profiles_HostTimeMatlab(end);
    else
        continue
    end
end
%leave NaNs as gaps so future data analysis doesn't include gaps as values
bcv(nanidx:length(bcv),:) = [];
if cmgidgaps(bcv) > 0
    [gaps,~,~,lgap] = cmgidgaps(bcv(:,1));
    disp(['There are ' num2str(gaps) ' gaps in ' fields{bci}])
    disp(['Gaps are ' num2str(lgap) ' samples long'])
end
DATA.BottomCheck_BottomDistance = bcv;
clearvars bcv data n m idx

%vecpro files also have a lot of metadata in the 'Data' structure. To
%account for these fields, I will copy these fields (that do not change)
%from the first file.
%first, find the fields not all ready concatenated
fnums = [nums bcti bci];
sfnums = setxor(1:length(fields),fnums);
disp('Writing ancillary data to structure')
for j = sfnums
    ind = find(strcmp(datanames(j),fields));
    DATA.(fields{ind}) = A(1).Data.(fields{ind});
    disp(['Copied ' fields{ind} ' to structure'])
end

%save out final structure
VPRO.Data = DATA;
VPRO.Config = A(1).Config;
clearvars DATA A time timev

%get file size info, provide a catch if final file is too large
info = whos('VPRO');
fsize = (info.bytes*1E-9)/4; %convert bytes to Gb, workspace variables are 4x actual file size
fprintf(['Structure VPRO is ' num2str(fsize) 'Gb\n'])

maxfilesize = 2; %Gb
if fsize > maxfilesize
    fprintf('WARNING: File size exceeds maximum alloted size of 2Gb\n')
    fprintf('Reduce output file size by reducing the number of input files\nExiting...\n')
    return
end

%save file
% save(fname,'VPRO','-v7.3')
% disp(['file saved as ' fname '.mat'])
% return

end
