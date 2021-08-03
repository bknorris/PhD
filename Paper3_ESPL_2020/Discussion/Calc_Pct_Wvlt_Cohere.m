%find % of wavelet files that are coherent vs. non-coherent
%GM_2018? BKN 2018
%
clear
%Define working paths
dpath = '\DataAnalysis\Paper3\';
ypath = {'Mekong_F2014';'Mekong_W2015'};
%load run file to tell program which files & VPs to load
rdir = 'g:\GradSchool\DataAnalysis\Paper3\ExperimentalDesign\';
fid = fopen([rdir 'AutoEventPCTrunfile.csv']);
rfile = textscan(fid,'%s%s%s%n%s%s','delimiter',',');
rdate = rfile{1};rstart = rfile{2};%dates corresponding to folder names; start time for crop
rstop = rfile{3};rvp = rfile{4};  %stop time for crop; vp #
wtdir = rfile{5};rarea = rfile{6};
start = datenum(strcat(rdate,{' '},rstart),'dd-mm-yy HH:MM:SS');
stop = datenum(strcat(rdate,{' '},rstop),'dd-mm-yy HH:MM:SS');
%initalize dataset
data.mud.fracevent = [];
data.fringe.fracevent = [];
data.forest.fracevent = [];
dfn = {'vpro1';'vpro2';'vpro3'};
%% Loop through folder structure
% for i = 1:2
% 	npath = ['f:\' ypath{i} dpath];
%     folders = dir([npath '\VPs\']);
%     folders = {folders(3:end).name};
%     for ii = 1:length(folders)
%         if any(strcmp(folders{ii},rdate))
%             vpid = strcmp(folders{ii},rdate);vpid = find(vpid);
%             disp(['Loading files from ' npath 'Wavelet\' folders{ii} '\'])
%             wfile  = dir([npath 'Wavelet\' folders{ii} '\','*.mat']); 
%         else
%             continue
%         end
%         for j = 1:length(vpid)
%             %% Load files
%             disp(['Processing ' dfn{rvp(vpid(j))}])
%             load([npath 'Wavelet\' folders{ii} '\' wfile(j).name]);
%             %% Filter by movement events
%             %find high coherence events (sig95 > 0.9) per bandwidth
%             threshold = 0.9;
%             [m,n] = size(wvlt.(wtdir{j}).sig95);
%             t = wvlt.(wtdir{j}).t;
%             events = zeros(m,n);
%             for k = 1:m
%                 events(k,:) = wvlt.(wtdir{j}).sig95(k,:) >= threshold;
%             end
%             %filter by coi
%             zid = zeros(m,n);
%             for k = 1:n
%                 zid(:,k) = wvlt.(wtdir{j}).period <= wvlt.x.coi(k);
%             end
%             %calculate fraction of wavelet analysis that is coherent and
%             %non coherent
%             events = events.*zid;
%             events(events==0) = NaN;
%             numcoi = nnz(zid);
%             isevent = length(find(~isnan(events)));
%             fracevent = isevent/numcoi;
%             %save data
%             if isempty(data.(rarea{vpid(j)}).fracevent)
%                 c = 1;
%             else
%                 c = length(data.(rarea{vpid(j)}).fracevent)+1;
%             end
%             data.(rarea{vpid(j)}).fracevent(c) = fracevent; 
%             clear wvlt
%         end
%     end
% end
load('F:\Mekong_W2015\DataAnalysis\Paper3\Wavelet\PCTcohere.mat')
fn = fieldnames(data);
for i = 1:3
    pct = data.(fn{i}).fracevent*100;
    CI = (std(pct')./sqrt(size(pct,2))).*1.96; %0.95 CI
    fprintf('%s percent coherent: %0.2f+/-%0.2f %%\n',fn{i},mean(pct),CI)
end