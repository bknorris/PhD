%This program automatically loads, pre-processes Vectrino data, then runs
%the cross-spectral wavelet analysis.
clear
%define working paths
dpath = '\DataAnalysis\Paper3\';
ypath = {'Mekong_F2014';'Mekong_W2015'};
%load run file to tell program which files & VPs to load
rdir = 'e:\MemoryStick\GradSchool\DataAnalysis\Paper3\ExperimentalDesign\';
fid = fopen([rdir 'AutoWLrunfile_v6.csv']);
rfile = textscan(fid,'%s%s%s%n%n%n','delimiter',',');
rdate = rfile{1};rstart = rfile{2};%dates corresponding to folder names; start time for crop
rstop = rfile{3};rvp = rfile{4};  %stop time for crop; vp #
rbin5 = rfile{5};rhead = rfile{6}; %use bins 1-5 for averaging; VP heading for cross & along-shore rot
start = datenum(strcat(rdate,{' '},rstart),'dd-mm-yy HH:MM:SS');
stop = datenum(strcat(rdate,{' '},rstop),'dd-mm-yy HH:MM:SS');
dfn = {'vpro1';'vpro2';'vpro3'};
for i = 1:2
    npath = ['e:\' ypath{i} dpath];
    folders = dir([npath '\VPs\']);
    folders = {folders(3:end).name};
    for ii = 1:length(folders)
        if any(strmatch(folders{ii},rdate))
            vpid = strmatch(folders{ii},rdate);
            disp(['Loading files from ' npath 'VPs\' folders{ii} '\'])
            disp(['Loading files from ' npath 'BottomTrack\' folders{ii} '\'])
            vfile = dir([npath 'VPs\' folders{ii} '\','*_Vels.mat']);
            bfile = dir([npath 'BottomTrack\' folders{ii} '\','*_bdtrace.mat']);
            bd = load([npath 'BottomTrack\' folders{ii} '\' bfile.name]);
        else
            continue
        end
        %just load the variables you need, saves memory
        for j = 1:length(vpid)
            disp(['Processing ' dfn{rvp(vpid(j))}])
            dat = load([npath 'VPs\' folders{ii} '\' vfile.name],dfn{rvp(vpid(j))});
            vpt = dat.(dfn{rvp(vpid(j))}).time;
            vpx = dat.(dfn{rvp(vpid(j))}).x;
            vpy = dat.(dfn{rvp(vpid(j))}).y;
            bdt = bd.(dfn{rvp(vpid(j))}).time;
            bds = bd.(dfn{rvp(vpid(j))}).bdist;
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
            %NOTE: some files have long data gaps. Files are cropped to
            %user-specified lengths in the 2-3 column of the run file (.csv).
            tid = find(vpt>=start(vpid(j))&vpt<=stop(vpid(j)));
            rxy = rxy(tid,:);
            vpt = vpt(tid,:);
            tid = find(bdt>=start(vpid(j))&bdt<=stop(vpid(j)));
            bds = bds(tid);
            bdt = bdt(tid);
            %Downsample velocity to 10Hz
            [vpt,idx]=unique(vpt);
            x10 = interp1(vpt,rxy(idx,1),bdt);nid = find(isnan(x10));
            disp(['Found ' num2str(length(nid)) ' NaNs in x10'])
            x10(nid) = 0;
            y10 = interp1(vpt,rxy(idx,2),bdt);nid = find(isnan(y10));
            disp(['Found ' num2str(length(nid)) ' NaNs in y10'])
            y10(nid) = 0;
            nid = find(isnan(bds));
            disp(['Found ' num2str(length(nid)) ' NaNs in bottom trace'])
            bds(nid) = 0;
            %Karin recommends detrending each t-s, Grinstead et al. 2004
            %recommends plotting the histograms of each t-s to see if they
            %are roughly gaussian. NOTE: I have done this, detrending seems
            %to make little difference for the velocities, but a big
            %difference (in terms of being closer to gaussian) for the
            %bottom trace t-s. 
            x10 = detrend(x10);y10 = detrend(y10);
            bds = detrend(bds);
            fs = 10;dt = 1/fs;
            time = (dt:dt:length(x10)*dt)/60; %time in minutes for x-axis plotting
            xt = [time' x10];
            yt = [time' y10];
            bt = [time' bds];
            [Rsq,period,scale,coi,sig95,Wxy,t,dt]=wtc(xt,bt,'S0',10/60);
            %Save data files for x-b
            wvlt.x.Rsq = Rsq;
            wvlt.x.period = period;
            wvlt.x.scale = scale;
            wvlt.x.coi = coi;
            wvlt.x.sig95 = sig95;
            wvlt.x.Wxy = Wxy;
            wvlt.x.t = t;
            wvlt.x.dt = dt;
            [Rsq,period,scale,coi,sig95,Wxy,t,dt]=wtc(yt,bt,'S0',10/60);
            %Save data files for y-b
            wvlt.y.Rsq = Rsq;
            wvlt.y.period = period;
            wvlt.y.scale = scale;
            wvlt.y.coi = coi;
            wvlt.y.sig95 = sig95;
            wvlt.y.Wxy = Wxy;
            wvlt.y.t = t;
            wvlt.y.dt = dt;
            save([npath 'Wavelet\' folders{ii} '\' dfn{rvp(vpid(j))} '_wvlt'],'wvlt','-v7.3')
            clear dat wvlt
        end
        clear bd
    end
end