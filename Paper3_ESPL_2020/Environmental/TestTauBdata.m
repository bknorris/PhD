%script compares old/new Taub data (as of 03/2019) to check for
%differences.
clear, close all
%define working paths
dpath = '\DataAnalysis\Paper3\';
ypath = {'Mekong_F2014';'Mekong_W2015'};
rdir = 'd:\GradSchool\DataAnalysis\Paper3\ExperimentalDesign\';
fid = fopen([rdir 'AutoTaucRunfile.csv']);
rfile = textscan(fid,'%s%s%s%n%n%s%n%n%n%n%n','delimiter',',');
rdate = rfile{1};%dates corresponding to folder names; start time for
rvp = rfile{4};  %stop time for crop; vp #
dfn = {'vpro1';'vpro2';'vpro3'};
dmin = zeros(3,1);
dmax = zeros(3,1);
for i = 1:2
    npath = ['e:\' ypath{i} dpath];
    folders = dir([npath '\VPs\']);
    folders = {folders(3:end).name};
    for ii = 1:length(folders)
        if any(strcmp(folders{ii},rdate))
            vpid = strcmp(folders{ii},rdate);
            vpid = find(vpid);
            figure
            %Just load the variables you need, saves memory
            for j = 1:length(vpid)
                vfile1 = dir([npath 'BedStress\' folders{ii} '\',...
                    '*' [dfn{rvp(vpid(j))} '.mat']]);
                vfile2 = dir([npath 'BedStress\' folders{ii} '\',...
                    '*' [dfn{rvp(vpid(j))} '_v2.mat']]);
                dat1 = load([vfile1.folder '\' vfile1.name]);
                dat2 = load([vfile2.folder '\' vfile2.name]);
                dmin(j) = min([nanmin(dat1.dat.tmax) nanmin(dat2.dat.tmax)]);
                dmax(j) = max([nanmax(dat1.dat.tmax) nanmax(dat2.dat.tmax)]);
                sp(j) = subplot(1,length(vpid),j);
                p(1) = plot(dat1.dat.time,dat1.dat.tmax,'b');hold on
                p(2) = plot(dat2.dat.time,dat2.dat.tmax,'r');
                legend(p,{'Original','Modified'})
                datetickzoom('x','dd HH:MM:SS','keepticks','keeplimits')
                xlabel('Time')
                ylabel('Shear Stress (N/m^2)')
                title(dfn{rvp(vpid(j))})
            end
            set(sp,'ylim',[min(dmin) max(dmax)])
            suplabel(folders{ii},'t');
        end
    end
end