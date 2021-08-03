%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Data processing script for the Nortek Vectrino II Profiler.
%Run this script in the folder that contains the raw data files.

%Script loads and concatenates VP files, then runs basic QC. Timebases are
%adjusted to local (VN) time. 

% This script was written by Benjamin K Norris, 2015
% University of Waikato, New Zealand

%Updated 16/05/2017
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear
close all

% toffset = [5 6 6 6 5 6 6 5]; %2014 time offset (hrs) used to adjust VP time to local time (different for each VP)
toffset = repmat(7,1,12); %2015 time offset (hrs) used to adjust VP time to local time (different for each VP)
year = '2015';

%find files
fdir = 'd:\Projects\Mekong_W2015\Data/Vectrino\';
tier1 = dir(fdir);
isub = [tier1(:).isdir];
tier1 = {tier1(isub).name}';
tier1 = tier1(3:end);
for i = 12 %:length(tier1)
    subfdir = [fdir tier1{i} '\'];
    tier2 = dir([subfdir '*.mat']);tier2 = {tier2.name};
    %now separate file list by VP (e.g. VP1, VP2, etc.)
    vps = {'VP1';'VP2';'VP3'};
    for ii = 1:3
        disp(['Processing files for ' tier1{i} ' ' vps{ii}])
        tic
        disp(['Run started at: ' datestr(now)])
        %%%
        files = strfind(upper(tier2),vps{ii});
        fid = find(~cellfun(@isempty,files));
        fileList = cell(1,length(fid));
        for id = 1:length(fid)
            fileList{id} = [subfdir tier2{fid(id)}];
        end
        %Concat VPs
        disp('Concatenating Vectrino Files')
        str = regexprep(tier1{i},'[^a-zA-Z0-9]','_');
        newFile = [str '_' vps{ii}];
        VPRO = concat_vecpro(fileList,newFile);
        
        %Run QC
        disp('Running QC processing')
        [VPRO.Data,VPRO.Config] = VPQC(VPRO.Data,VPRO.Config,1,1,0);        %run QC on data
        
        %Adjust HostTimeMatlab to local VN time
        if strcmp(year,'2014')
            tc = -1*datenum(0,0,0,toffset(i),0,0);
        elseif strcmp(year,'2015')
            tc = datenum(0,0,0,toffset(i),0,0);
        end
        VPRO.Data.Time = VPRO.Data.Time+tc;
        VPRO.Data.Profiles_HostTimeMatlab = VPRO.Data.Profiles_HostTimeMatlab+tc;
        VPRO.Data.BottomCheck_HostTimeMatlab = VPRO.Data.BottomCheck_HostTimeMatlab+tc;
        
        save([subfdir newFile],'VPRO','-v7.3')
        disp(['Basic Quality Controls completed in: ' num2str(toc/60) ' minutes'])
        clear VPRO
    end
end
    