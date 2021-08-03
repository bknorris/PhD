%Check that trends shown in events for the BdPDF_alldata Probability
%Density Functions are consistent across all of the bed level data.
%
clear
%% Load and calculate bed level statistics for all events
%define working paths
dpath = '\DataAnalysis\Paper3\';
ypath = 'Mekong_W2015';
%load run file to tell program which files & VPs to load
rdir = 'e:\MemoryStick\GradSchool\DataAnalysis\Paper3\ExperimentalDesign\';
fname = 'AutoMLRrunfile_ebb.csv';
fid = fopen([rdir fname]);
rfile = textscan(fid,'%s%n%s%n%n%n','delimiter',',');
rdate = rfile{1};
rexn = rfile{2};
rexp = rfile{3};
rvph = rfile{4};
rd = rfile{5};
rD = rfile{6};
dat = struct();
fields = {'deltbd';'minbd';'maxbd';'medbd';'netbd'};
flds = flipud(unique(rexp));
for i = 1:length(flds)
    for ii = 1:length(fields)
        dat.(flds{i}).(fields{ii}) = [];
    end
end
npath = ['e:\' ypath dpath];
folders = dir([npath '\BottomTrack\']);
folders = {folders(3:end).name};
for ii = 1:length(folders)
    if any(strcmp(folders{ii},rdate))
        disp(['Loading files from ' npath 'BottomTrack\' folders{ii} '\'])
        vpid = strcmp(folders{ii},rdate);vpid = find(vpid);
        file = dir([npath 'BottomTrack\' folders{ii} '\','*bdtrace.mat']);
        data = load([npath 'BottomTrack\' folders{ii} '\',file.name]);
    else
        continue
    end
    dfn = fieldnames(data);
    for j = 1:length(dfn)
        exp = rexp{vpid(j)};
        bdist = data.(dfn{j}).bdist;
        bdstart = find(~isnan(data.(dfn{j}).bdist),1,'first');
        bdend = find(~isnan(data.(dfn{j}).bdist),1,'last');
        bdist = bdist(bdstart)-bdist;
        %calculate elevation stats for individual time series
        deltbd = (bdist(end)-bdist(1))/length(bdist);
        minbd = nanmin(bdist);
        maxbd = nanmax(bdist);
        medbd = nanmedian(bdist);
        netbd = bdist(bdend)-bdist(bdstart);
        dat.(exp).deltbd = [dat.(exp).deltbd deltbd];
        dat.(exp).minbd = [dat.(exp).minbd minbd];
        dat.(exp).maxbd = [dat.(exp).maxbd maxbd];
        dat.(exp).medbd = [dat.(exp).medbd medbd];
        dat.(exp).netbd = [dat.(exp).netbd netbd];
    end
end
%also load data from 09/03/15
load('e:\Mekong_W2015\DataAnalysis\Paper3\BottomTrack\BDtrace_floodebb.mat')
dfn = fieldnames(data);
%just grab data from VPRO2
%crop time based on flood vs. ebbexp = 'fringe';
time = data.(dfn{2}).time;
bdist = data.(dfn{2}).ble;
hightide = datenum(2015,03,09,04,40,00);
if contains(fname,'_flood')
    tidx = find(time >= time(1) & time <= hightide);
elseif contains(fname, '_ebb')
    tidx = find(time >= hightide & time <= time(end));
end
bdist = bdist(tidx);
%calculate elevation stats for individual time series
deltbd = (bdist(end)-bdist(1))/length(bdist);
minbd = nanmin(bdist);
maxbd = nanmax(bdist);
medbd = nanmedian(bdist);
netbd = bdist(end)-bdist(find(~isnan(bdist),1,'first'));
dat.fringe.deltbd = [dat.fringe.deltbd deltbd];
dat.fringe.minbd = [dat.fringe.minbd minbd];
dat.fringe.maxbd = [dat.fringe.maxbd maxbd];
dat.fringe.medbd = [dat.fringe.medbd medbd];
dat.fringe.netbd = [dat.fringe.netbd netbd];
        
dfn = fieldnames(dat);
for i = 1:length(dfn)
    fprintf('Min BLE across all events for %s: %0.1f mm\n',dfn{i},min(dat.(dfn{i}).minbd)*1000)
    fprintf('Median BLE across all events for %s: %0.1f mm\n',dfn{i},median(dat.(dfn{i}).medbd)*1000)
    fprintf('Max BLE across all events for %s: %0.1f mm\n',dfn{i},max(dat.(dfn{i}).maxbd)*1000)
    fprintf('Net BLE across all events for %s: %0.1f mm\n\n',dfn{i},nanmean(dat.(dfn{i}).netbd)*1000)
end

