%Plot turbulence and bed trace to look for instances where turb
%measurements don't exist. Will need to crop time series at this point.

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
        tfile = dir([tpath,'*_TKE.mat']);
        bfile = dir([bpath,'*bdtrace.mat']);
        tke = load([tpath tfile.name]);
        bd = load([bpath bfile.name]);
        fn = fieldnames(tke);
        f(ii) = figure(ii);
        maxlim = max(rhab(inid));
        for j = 1:length(fn)
            bh = rhab(inid(j))-bd.(fn{j}).bdist;
            bdt = bd.(fn{j}).time;
            tdt = tke.(fn{j}).time;
            vh = rhab(inid(j))-0.04-linspace(0.001,0.03,30);
            [~,bid] = unique(bdt);
            bhadj = interp1(bdt(bid),bh(bid),tdt)+0.005; %buffer zone of bed (5 mm)
            
            %Plot routine
            sp(j) = subplot(1,3,j);
            imagesc(tdt,vh,(tke.(fn{j}).z1.E+tke.(fn{j}).z2.E)./2)
            hold on
            plot(tdt,bhadj,'r')
            datetickzoom('x','HH:MM')
            set(gca,'xlim',[min(tdt) max(tdt)],'ylim',[0 maxlim],'ydir','normal')
            title([rdate{inid(j)} ' ' fn{j}])
        end
        cb = colorbar;
        set(sp(1),'position',[0.12 0.12 0.2 0.8])
        set(sp(2),'position',[0.38 0.12 0.2 0.8])
        set(sp(3),'position',[0.65 0.12 0.2 0.8])
        set(cb,'position',[0.88 0.12 0.015 0.8])
        ylabel(sp(1),'Height Above Bottom (m)')
        xlabel(sp(2),'Time of Day (HH:MM)')
        ylabel(cb,'TKE Dissipation Rate (W/m^2)')
        prettyfigures('text',8,'labels',10,'box',1)
        sfpath = 'e:\Mekong_W2015\DataAnalysis\Paper3\Turbulence\Figures\';
%         export_fig([sfpath rdate{inid(1)} '_TKEts'])
        close
    end
end