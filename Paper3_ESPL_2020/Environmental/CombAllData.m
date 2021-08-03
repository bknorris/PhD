clear,close all
%define working paths
dpath = '\DataAnalysis\Paper3\';
ypath = {'Mekong_F2014';'Mekong_W2015'};
%load run file to tell program which files & VPs to load
rdir = 'd:\MemoryStick\GradSchool\DataAnalysis\Paper3\ExperimentalDesign\';
fid = fopen([rdir 'AutoMLRrunfile_ebb.csv']);
rfile = textscan(fid,'%s%n%s%n%n%n','delimiter',',');
rdate = rfile{1};
rexn = rfile{2};
rexp = rfile{3};
rvph = rfile{4};
rd = rfile{5};
rD = rfile{6};
%initialize structure
dat = struct();
fields = {'time';'umag';'uvmag';'usqd';'bdist';...
    'depth';'Hs';'ubr';'tmax';'eps'};
flds = flipud(unique(rexp));
for i = 1:length(flds)
    for ii = 1:length(fields)
        dat.(flds{i}).(fields{ii}) = [];
    end
    dat.(flds{i}).KC = [];
    dat.(flds{i}).d = [];
    dat.(flds{i}).D = [];
    dat.(flds{i}).expnum = [];
end
for i = 1:2
    npath = ['d:\' ypath{i} dpath];
    folders = dir([npath '\CmbData\']);
    folders = {folders(3:end).name};
    for ii = 1:length(folders)
        if any(strcmp(folders{ii},rdate))
            vpid = strcmp(folders{ii},rdate);vpid = find(vpid);
            disp(['Loading files from ' npath 'CmbData\' folders{ii} '\'])
            load([npath 'CmbData\' folders{ii} '\','AllEventData_v3.mat']);
        else
            continue
        end
        dfn = fieldnames(data);
        for j = 1:length(dfn)
            ffn = fieldnames(data.(dfn{j}));
            [~,ib] = intersect(ffn,fields);ib = sort(ib);
            exp = rexp{vpid(j)};
            for k = 1:length(ib)
                tmp = data.(dfn{j}).(ffn{ib(k)});
                tmp = tmp(~isnan(tmp));
                if strcmp(ffn{ib(k)},'bdist')
                    bst = find(~isnan(data.(dfn{j}).bdist),1,'first');
                    tmp = tmp-rvph(vpid(j));
                end
                dat.(exp).(fields{k}) = [dat.(exp).(fields{k}) ...
                    tmp];
            end
            %Calculate KC numbers
            ubr = data.(dfn{j}).ubr;
            tr = data.(dfn{j}).Tr;
            if rd(vpid(j)) == 0
                KCd = NaN(1,length(tr));
            else
                KCd = (ubr.*tr)/rd(vpid(j));
            end
%             KCD = (ubr.*tr)/rD(vpid(j));
            dat.(exp).KC = [dat.(exp).KC KCd];
            dat.(exp).d = [dat.(exp).d rd(vpid(j))];
            dat.(exp).D = [dat.(exp).D rD(vpid(j))];
            dat.(exp).expnum = [dat.(exp).expnum repmat(rexn(vpid(j)),1,length(tr))];
        end   
    end
end

sdir = 'd:\Mekong_W2015\DataAnalysis\Paper3\CmbData\';
save([sdir 'CombAllData_ebb_v3'],'dat','-v7.3')