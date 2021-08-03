clear,close all
%define working paths
dpath = '\DataAnalysis\Paper3\';
ypath = {'Mekong_F2014';'Mekong_W2015'};
%load run file to tell program which files & VPs to load
rdir = 'e:\MemoryStick\GradSchool\DataAnalysis\Paper3\ExperimentalDesign\';
fid = fopen([rdir 'AutoMLRrunfile_flood.csv']);
rfile = textscan(fid,'%s%n%s%n%n%n','delimiter',',');
rdate = rfile{1};
rexn = rfile{2};
rexp = rfile{3};
%initialize structure
dat = struct();
fields = {'eventl';'deltbd';'bdm20';'bdmnm';'usqd';...
    'phase';'sigh';'depth';'orbwv';'taub';'eps20';'epsnm'};
flds = flipud(unique(rexp));
for i = 1:length(flds)
    for ii = 1:length(fields)
        dat.(flds{i}).wave.(fields{ii}) = [];
        dat.(flds{i}).ig.(fields{ii}) = [];
    end
end
for i = 1:2
    npath = ['e:\' ypath{i} dpath];
    folders = dir([npath '\VPs\']);
    folders = {folders(3:end).name};
    for ii = 1:length(folders)
        if any(strcmp(folders{ii},rdate))
            vpid = strcmp(folders{ii},rdate);vpid = find(vpid);
            disp(['Loading files from ' npath 'Wavelet\' folders{ii} '\'])
            load([npath 'Wavelet\' folders{ii} '\','waveIGevents_v5.mat']);
        else
            continue
        end
        dfn = fieldnames(eventd);
        for j = 1:length(dfn)
            %WAVE
            ffn = fieldnames(eventd.(dfn{j}).wave);
            [~,ib] = intersect(ffn,fields);ib = sort(ib);
            exp = rexp{vpid(j)};
%             disp(['Wave ' exp ' ' dfn{j}])
            for k = 1:length(ib)
                tmp = eventd.(dfn{j}).wave.(ffn{ib(k)});
                tmp = tmp(~isnan(tmp));
                if strcmp(ffn{ib(k)},'bdmed')
                    tmp = tmp-eventd.(dfn{j}).wave.bdst;
                end
                dat.(exp).wave.(fields{k}) = [dat.(exp).wave.(fields{k});...
                    tmp];
            end
            tmp = repmat(rexn(vpid(j)),length(tmp),1);
            %IG
            ffn = fieldnames(eventd.(dfn{j}).ig);
            [~,ib] = intersect(ffn,fields);ib = sort(ib);
            exp = rexp{vpid(j)};
%             disp(['IG ' exp ' ' dfn{j}])
            for k = 1:length(ib)
                tmp = eventd.(dfn{j}).ig.(ffn{ib(k)});
                tmp = tmp(~isnan(tmp));
                if strcmp(ffn{ib(k)},'bdmed')
                    tmp = tmp-eventd.(dfn{j}).ig.bdst;
                end
                dat.(exp).ig.(fields{k}) = [dat.(exp).ig.(fields{k});...
                    tmp];
            end
            tmp = repmat(rexn(vpid(j)),length(tmp),1);
        end
    end
end

sdir = 'e:\Mekong_W2015\DataAnalysis\Paper3\Wavelet\';
save([sdir 'CombEventData_flood_v5'],'dat','-v7.3')