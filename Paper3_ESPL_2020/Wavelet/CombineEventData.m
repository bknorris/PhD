clear,close all
%define working paths
dpath = '\DataAnalysis\Paper3\';
ypath = {'Mekong_F2014';'Mekong_W2015'};
%load run file to tell program which files & VPs to load
rdir = 'f:\GradSchool\DataAnalysis\Paper3\ExperimentalDesign\';
fid = fopen([rdir 'AutoMLRrunfile.csv']);
rfile = textscan(fid,'%s%s','delimiter',',');
rdate = rfile{1};rexp = rfile{2};
%Define model parameters
eventl = NaN(11000,length(rdate));
bdvar = NaN(11000,length(rdate));
deltbd = NaN(11000,length(rdate));
bdmed = NaN(11000,length(rdate));
Hs = NaN(11000,length(rdate));
H = NaN(11000,length(rdate));
uorb = NaN(11000,length(rdate));
umag = NaN(11000,length(rdate));
for i = 1:2
    npath = ['g:\' ypath{i} dpath];
    folders = dir([npath '\VPs\']);
    folders = {folders(3:end).name};
    for ii = 1:length(folders)
        if any(strcmp(folders{ii},rdate))
            vpid = strcmp(folders{ii},rdate);vpid = find(vpid);
            disp(['Loading files from ' npath 'Wavelet\' folders{ii} '\'])
            load([npath 'Wavelet\' folders{ii} '\','waveIGevents.mat']);
        else
            continue
        end
        dfn = fieldnames(eventd);
        for j = 1:length(dfn)
            %to start, I'll combine wave and IG events. In later renditions
            %of this script I may separate them for individual analysis.
            %%%
            tmp = [eventd.(dfn{j}).wave.eventl; eventd.(dfn{j}).ig.eventl];
            tmp = tmp(~isnan(tmp));
            eventl(1:length(tmp),vpid(j)) = tmp;
            tmp = [eventd.(dfn{j}).wave.bdvar; eventd.(dfn{j}).ig.bdvar];
            tmp = tmp(~isnan(tmp));
            bdvar(1:length(tmp),vpid(j)) = tmp;
            tmp = [eventd.(dfn{j}).wave.deltbd; eventd.(dfn{j}).ig.deltbd];
            tmp = tmp(~isnan(tmp));
            deltbd(1:length(tmp),vpid(j)) = tmp;
            tmp = [eventd.(dfn{j}).wave.bdmed; eventd.(dfn{j}).ig.bdmed];
            tmp = tmp(~isnan(tmp))-eventd.(dfn{j}).wave.bdst;
            bdmed(1:length(tmp),vpid(j)) = tmp;
            tmp = [eventd.(dfn{j}).wave.sigh; eventd.(dfn{j}).ig.sigh];
            tmp = tmp(~isnan(tmp));
            Hs(1:length(tmp),vpid(j)) = tmp;
            tmp = [eventd.(dfn{j}).wave.depth; eventd.(dfn{j}).ig.depth];
            tmp = tmp(~isnan(tmp));
            H(1:length(tmp),vpid(j)) = tmp;
            tmp = [eventd.(dfn{j}).wave.orbwv; eventd.(dfn{j}).ig.orbwv];
            tmp = tmp(~isnan(tmp));
            uorb(1:length(tmp),vpid(j)) = tmp;
            tmp = [eventd.(dfn{j}).wave.umag; eventd.(dfn{j}).ig.umag];
            tmp = tmp(~isnan(tmp));
            umag(1:length(tmp),vpid(j)) = tmp;
        end
    end
end
%initialize structure
dat = struct();
fields = {'eventl';'bdvar';'deltbd';'bdmed';'Hs';'H';'uorb';'umag'};
flds = flipud(unique(rexp));
for i = 1:length(flds)
    for ii = 1:length(fields)
        dat.(flds{i}).(fields{ii}) = NaN(11000,length(rdate));
    end
end
%fill structure with data in correct position
for i = 1:length(rdate)
    exp = rexp{i}; %position of experiment (mudflat/fringe/forest)
    dat.(exp).eventl(:,i) = eventl(:,i);
    dat.(exp).bdvar(:,i) = bdvar(:,i);
    dat.(exp).deltbd(:,i) = deltbd(:,i);
    dat.(exp).bdmed(:,i) = bdmed(:,i);
    dat.(exp).Hs(:,i) = Hs(:,i);
    dat.(exp).H(:,i) = H(:,i);
    dat.(exp).uorb(:,i) = uorb(:,i);
    dat.(exp).umag(:,i) = umag(:,i);
end
%remove column NaNs
for i = 1:length(flds)
    for ii = 1:length(fields)
        tmp = dat.(flds{i}).(fields{ii});
        tmp(~any(~isnan(tmp),2),:) = []; %trailing nans
        tmp(:,~any(~isnan(tmp))) = []; %col nans
        [m,n] = size(tmp);
        tmp = reshape(tmp,m*n,1);
        dat.(flds{i}).(fields{ii}) = tmp;
    end
end
sdir = 'g:\Mekong_W2015\DataAnalysis\Paper3\Wavelet\';
save([sdir 'CombEventData'],'dat','-v7.3')