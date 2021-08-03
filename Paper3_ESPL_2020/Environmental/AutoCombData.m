%This script is a copy of EventDuration.m, and simply loads and
%concatenates the individual experiment files into groupings, in order to
%make a complete dataset of 
clear
%Define working paths
dpath = '\DataAnalysis\Paper3\';
ypath = {'Mekong_F2014';'Mekong_W2015'};
%load run file to tell program which files & VPs to load
rdir = 'd:\MemoryStick\GradSchool\DataAnalysis\Paper3\ExperimentalDesign\';
fid = fopen([rdir 'AutoEventDrunfile.csv']);
rfile = textscan(fid,'%s%s%s%n%n%n%s%s','delimiter',',');
rdate = rfile{1};rstart = rfile{2};%dates corresponding to folder names; start time for crop
rstop = rfile{3};rvp = rfile{4};  %stop time for crop; vp #
rbin5 = rfile{5};rhead = rfile{6}; %use bins 1-5 for averaging; VP heading for cross & along-shore rot
wtdir = rfile{7};aqfil = rfile{8}; %wavelet direction ('x' or 'y'), aquadopp to load
start = datenum(strcat(rdate,{' '},rstart),'dd-mm-yy HH:MM:SS');
stop = datenum(strcat(rdate,{' '},rstop),'dd-mm-yy HH:MM:SS');
dfn = {'vpro1';'vpro2';'vpro3'};
tic
%% Loop through folder structure
for i = 1:2
	npath = ['d:\' ypath{i} dpath];
    folders = dir([npath '\VPs\']);
    folders = {folders(3:end).name};
    for ii = 1:length(folders)
        if any(strcmp(folders{ii},rdate))
            vpid = strcmp(folders{ii},rdate);vpid = find(vpid);
            disp(['Loading files from ' npath 'VPs\' folders{ii} '\'])
            disp(['Loading files from ' npath 'BottomTrack\' folders{ii} '\'])
            disp(['Loading files from ' npath 'Wavelet\' folders{ii} '\'])
            disp(['Loading files from ' npath 'WaveStats\' folders{ii} '\'])
            disp(['Loading files from ' npath 'BedStress\' folders{ii} '\'])
            disp(['Loading files from ' npath 'Turbulence\' folders{ii} '\'])
            vfile = dir([npath 'VPs\' folders{ii} '\','*_Vels.mat']);
            wfile  = dir([npath 'Wavelet\' folders{ii} '\','*.mat']); 
            bfile = dir([npath 'BottomTrack\' folders{ii} '\','*_bdtrace.mat']);
            afile = dir([npath 'WaveStats\' folders{ii} '\','*_v2.mat']);afile = {afile.name};
            tfile = dir([npath 'Turbulence\' folders{ii} '\','*.mat']);
            sfile = dir([npath 'BedStress\' folders{ii} '\','*_v2.mat']);
            bd = load([npath 'BottomTrack\' folders{ii} '\' bfile.name]);
        else
            continue
        end
        data = struct();
        for j = 1:length(vpid)
            %% Load files
            disp(['Processing ' dfn{rvp(vpid(j))}])
            load([npath 'Wavelet\' folders{ii} '\' wfile(j).name]);
            aqid = strfind(afile, aqfil{vpid(j)});aqid = find(~cellfun(@isempty,aqid));
            load([npath 'WaveStats\' folders{ii} '\' afile{aqid}])
            tke = load([npath 'Turbulence\' folders{ii} '\' tfile.name],dfn{rvp(vpid(j))});
            tke = tke.(dfn{rvp(vpid(j))});
            load([npath 'BedStress\' folders{ii} '\' sfile(j).name]);     
            tau = dat;clear dat
            dat = load([npath 'VPs\' folders{ii} '\' vfile.name],dfn{rvp(vpid(j))});
            %Define Variables
            vpt = dat.(dfn{rvp(vpid(j))}).time;
            vpx = dat.(dfn{rvp(vpid(j))}).x;
            vpy = dat.(dfn{rvp(vpid(j))}).y;
            bdt = bd.(dfn{rvp(vpid(j))}).time;
            bds = bd.(dfn{rvp(vpid(j))}).bdist;
            Hs = wave.Hs;
            H = wave.h;
            tr = wave.Tr;
            ubr = wave.ubr;
            wvt = wave.time2;
            ept = tke.time;
            eps = nanmean((tke.z1.E+tke.z2.E)./2);
            eps(eps>1E-3) = 0; %remove anomalous large epsilon estimates
            eps(isnan(eps)) = 0;
            %Some of the tmax estimates are corrupted by measurements too
            %close to the bed, set these equal to tauw
            tmax = tau.tmax; %should be length(wvt)
            if nnz(isnan(tmax)) > round(length(tmax)/2)
                tmax = tau.tauw;
            end
            tmax(isnan(tmax)) = 0;
            %% Rotate VP velocities to cross & along shore
            if rbin5(vpid(j)) == 1 %if VP is close to ground, use bins 1-5
                bins = 1:5;
            else
                bins = 13:17;
            end
            heading = rhead(vpid(j));th = heading*pi/180;
            R = [cos(th) -sin(th); sin(th) cos(th)];
            vpx = nanmean(vpx(:,bins),2);
            vpy = nanmean(vpy(:,bins),2);
            vpxy = [vpx vpy];rxy = zeros(size(vpxy));
            for jj = 1:length(vpxy)
                rxy(jj,:) = vpxy(jj,:)*R;end %rotate x, y to cross & along-shore
            %% (Down)Upsample data to 10 Hz
            %NOTE: some files have long data gaps. Files are cropped to
            %user-specified lengths in the 2-3 column of the run file (.csv).
            tid = find(vpt>=start(vpid(j))&vpt<=stop(vpid(j)));
            rxy = rxy(tid,:);
            vpt = vpt(tid,:);
            tid = find(bdt>=start(vpid(j))&bdt<=stop(vpid(j)));
            bds = bds(tid);
            bdt = bdt(tid);
            %Downsample velocity to 1Hz
            [vpt,idx]=unique(vpt);
            [ept,idx2]=unique(ept);
            x1 = interp1(vpt,rxy(idx,1),ept);nid = find(isnan(x1));
            disp(['Found ' num2str(length(nid)) ' NaNs in x1'])
            x1(nid) = [];
            y1 = interp1(vpt,rxy(idx,2),ept);nid = find(isnan(y1));
            disp(['Found ' num2str(length(nid)) ' NaNs in y1'])
            y1(nid) = [];
            ept(nid) = [];eps(nid) = [];
            %Interpolate wave parameters to 1Hz
            H1 = interp1(wvt,H,ept);
            Hs1 = interp1(wvt,Hs,ept);
            Tr1 = interp1(wvt,tr,ept);
            Om1 = interp1(wvt,ubr,ept);
            %Interpolate stress & bed level parameters to 1Hz
            tb1 = interp1(wvt,tmax,ept);
            [bdt,idx] = unique(bdt);
            bs1 = interp1(bdt,bds(idx),ept);
            %% Save data to structure
            data.(dfn{rvp(vpid(j))}).time = ept;
            data.(dfn{rvp(vpid(j))}).x = x1;
            data.(dfn{rvp(vpid(j))}).y = y1;
            data.(dfn{rvp(vpid(j))}).umag = sqrt(x1.^2);
            data.(dfn{rvp(vpid(j))}).uvmag = sqrt(x1.^2+y1.^2);
            data.(dfn{rvp(vpid(j))}).usqd = x1.^2;
            data.(dfn{rvp(vpid(j))}).ucub = x1.^3;
            data.(dfn{rvp(vpid(j))}).bdist = bs1;
            data.(dfn{rvp(vpid(j))}).depth = H1;
            data.(dfn{rvp(vpid(j))}).Hs = Hs1;
            data.(dfn{rvp(vpid(j))}).Tr = Tr1;
            data.(dfn{rvp(vpid(j))}).ubr = Om1;
            data.(dfn{rvp(vpid(j))}).tmax = tb1;
            data.(dfn{rvp(vpid(j))}).eps = eps;
            clear dat wvlt tke tau
            save([npath 'CmbData\' folders{ii} '\' 'AllEventData_v3'],'data','-v7.3')
        end
    end
end
fprintf('Data loaded, processed and concatenated in %0.2f minutes\n',toc/60)