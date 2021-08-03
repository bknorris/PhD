%script loads new Taub data (as of 03/2019) to see how many near-bed points
%were fit across all experiments.

clear, close all
%define working paths
dpath = '\DataAnalysis\Paper3\';
ypath = {'Mekong_F2014';'Mekong_W2015'};
rdir = 'd:\GradSchool\DataAnalysis\Paper3\ExperimentalDesign\';
fid = fopen([rdir 'AutoTaucRunfile.csv']);
rfile = textscan(fid,'%s%s%s%n%n%s%n%n%n%n%n','delimiter',',');
rdate = rfile{1};rstart = rfile{2};%dates corresponding to folder names; start time for crop
rstop = rfile{3};rvp = rfile{4};  %stop time for crop; vp #
rhead = rfile{5};raqdp = rfile{6}; %VP heading for cross & along-shore rot;aquadopp file name to load
rsmpr= rfile{7};rhab = rfile{8}; %aqdp sample rate (Hz),vp height above bed (m)
rd50 = rfile{9};rvegd = rfile{10}; %d50, >0 if veg, 0 if no veg
rfd50 = rfile{11}; %force use of d50 for calculation of z0
dfn = {'vpro1';'vpro2';'vpro3'};
rsqmin = NaN(length(rdate),1);
rsqmax = NaN(length(rdate),1);
count = 0;
for i = 1:2
    npath = ['e:\' ypath{i} dpath];
    folders = dir([npath '\VPs\']);
    folders = {folders(3:end).name};
    for ii = 1:length(folders)
        if any(strcmp(folders{ii},rdate))
            vpid = strcmp(folders{ii},rdate);
            vpid = find(vpid);
            %Just load the variables you need, saves memory
            for j = 1:length(vpid)
                count = count+1;
                vfile = dir([npath 'BedStress\' folders{ii} '\',...
                    '*' [dfn{rvp(vpid(j))} '_v2.mat']]);
                load([vfile.folder '\' vfile.name]);
                rsqmin(count) = nanmin(dat.rsq)+5; %5 samples was 'minsamp' in AutoTaucVPs.m
                rsqmax(count) = nanmax(dat.rsq)+5;
            end
        end
    end
end
fprintf('Minimum number of samples utilized: %0.0f\n',nanmin(rsqmin))
fprintf('Maximum number of samples utilized: %0.0f\n',nanmax(rsqmax))