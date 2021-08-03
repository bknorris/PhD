% This program automatically loads, and crops the bed out of the turbulence
% estimates
%
%
% BKN, UoW 2017
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear
%define working paths
dpath = '\DataAnalysis\Paper3\';
ypath = {'Mekong_F2014';'Mekong_W2015'};
%load run file to tell program which files & VPs to load
rdir = 'e:\MemoryStick\GradSchool\DataAnalysis\Paper3\ExperimentalDesign\';
fid = fopen([rdir 'AutoTKErunfile_v2.csv']);
rfile = textscan(fid,'%s%s%s%s%n%n','delimiter',',');
rexpt = rfile{1};rdate = rfile{2}; %experiment name, experiment date
rstart = rfile{3};rstop = rfile{4}; %start time, stop time (based on VPs)
rinst = rfile{5};rhab = rfile{6}; %instrument name in folder, instrument height above bed
start = datenum(strcat(rdate,{' '},rstart),'dd-mm-yy HH:MM:SS');
stop = datenum(strcat(rdate,{' '},rstop),'dd-mm-yy HH:MM:SS');
stop = stop+datenum(0,0,0,0,3,0); %account for windowed indexing
for i = 1:2
    npath = ['e:\' ypath{i} dpath 'BottomTrack\'];
    folders = dir(npath);
    folders = {folders(3:end).name};
    for ii = 1:length(folders)
        if any(strcmp(folders{ii},rdate))
            inid = strcmp(folders{ii},rdate);
            inid = find(inid);
            tpath = ['e:\' ypath{i} dpath 'Turbulence\' folders{ii} '\'];
            bpath = [npath folders{ii} '\'];
            disp(['Loading files from ' tpath])
            disp(['Loading files from ' bpath])
        else
            continue
        end
        tfile = dir([tpath,'*_TKE_v2.mat']);
        bfile = dir([bpath,'*bdtrace.mat']);
        tke = load([tpath tfile.name]);
        bd = load([bpath bfile.name]);
        fn = fieldnames(tke);
        for j = 1:length(fn)
            bh = rhab(inid(j))-bd.(fn{j}).bdist;
            bdt = bd.(fn{j}).time;
            tdt = tke.(fn{j}).time;
            vh = rhab(inid(j))-0.04-linspace(0.001,0.03,30);
            [~,bid] = unique(bdt);
            bhadj = interp1(bdt(bid),bh(bid),tdt)+0.005; %buffer zone of bed (5 mm)
            for jj = 1:length(tdt)
                id = find(bhadj(jj) >= vh,1,'first');
                if isempty(id)
                    continue
                else
                    tke.(fn{j}).z1.E(id:end,jj) = 0;
                    tke.(fn{j}).z1.N(id:end,jj) = 0;
                    tke.(fn{j}).z2.E(id:end,jj) = 0;
                    tke.(fn{j}).z2.N(id:end,jj) = 0;
                end
            end
        end
        sfname = [rexpt{inid(j)} '_' rdate{inid(j)}(1:2) '_TKE'];
        disp(['Saving ' tpath sfname '_v2.mat'])
        save([tpath sfname],'-struct','tke','-v7.3')
    end
end
