%Try to extract the event duration from individual frequency bands of the
%wavelet analysis (i.e., how long a high-coherence event lasts). Plot this
%against the normalized bed level elevation variance (e.g., Staudt et al.
%2017)
%10-11-17: I have added code to load AQDP files and crop by event length.
%04-12-17: I have added dissipation rate and bed shear stress
%
%
clear
%Define working paths
dpath = '\DataAnalysis\Paper3\';
ypath = {'Mekong_F2014';'Mekong_W2015'};
%load run file to tell program which files & VPs to load
rdir = 'e:\GradSchool\DataAnalysis\Paper3\ExperimentalDesign\';
fid = fopen([rdir 'AutoEventDrunfile.csv']);
rfile = textscan(fid,'%s%s%s%n%n%n%s%s','delimiter',',');
rdate = rfile{1};rstart = rfile{2};%dates corresponding to folder names; start time for crop
rstop = rfile{3};rvp = rfile{4};  %stop time for crop; vp #
rbin5 = rfile{5};rhead = rfile{6}; %use bins 1-5 for averaging; VP heading for cross & along-shore rot
wtdir = rfile{7};aqfil = rfile{8}; %wavelet direction ('x' or 'y'), aquadopp to load
start = datenum(strcat(rdate,{' '},rstart),'dd-mm-yy HH:MM:SS');
stop = datenum(strcat(rdate,{' '},rstop),'dd-mm-yy HH:MM:SS');
dfn = {'vpro1';'vpro2';'vpro3'};
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
            afile = dir([npath 'WaveStats\' folders{ii} '\','*.mat']);afile = {afile.name};
            tfile = dir([npath 'Turbulence\' folders{ii} '\','*.mat']);
            sfile = dir([npath 'BedStress\' folders{ii} '\','*.mat']);
            bd = load([npath 'BottomTrack\' folders{ii} '\' bfile.name]);
        else
            continue
        end
        t1 = tic;
        eventd = struct();
        for j = 1:length(vpid)
            %% Load files
            t2 = tic;
            disp(['Processing ' dfn{rvp(vpid(j))}])
            disp('Loading data...')
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
            disp('Rotating velocities to across and along-shore...')
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
            disp('Decimating velocity to 10 Hz')
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
            %Interpolate wave parameters to 10Hz
            H10 = interp1(wvt,H,bdt);
            Hs10 = interp1(wvt,Hs,bdt);
            Om10 = interp1(wvt,ubr,bdt);
            %Interpolate stress & turbulence parameters to 10Hz
            tb10 = interp1(wvt,tmax,bdt);
            [~,idx] = unique(ept);
            ep10 = interp1(ept(idx),eps(idx),bdt,'linear','extrap');
            fprintf('Basic data processing completed in %0.2f minutes\n',toc(t2)/60)
            %% Filter by movement events
            %find high coherence events (sig95 > 0.9) per bandwidth
            threshold = 0.9;
            [m,n] = size(wvlt.(wtdir{j}).sig95);
            t = wvlt.x.t;
            events = zeros(m,n);
            for k = 1:m
                events(k,:) = wvlt.(wtdir{j}).sig95(k,:) >= threshold;
            end
            %filter by coi
            zid = zeros(m,n);
            for k = 1:n
                zid(:,k) = wvlt.(wtdir{j}).period <= wvlt.(wtdir{j}).coi(k);
            end
            events = events.*zid;
            %% do some 'image' processing; WARNING: Requires image processing toolbox!!!
            eventGroups = bwlabel(events,8);
            maxNumEvnt = max(max(eventGroups));
            %initialize variables to save out
            nevent = NaN(maxNumEvnt,1); %# timesteps in event
            deltabd = NaN(maxNumEvnt,1); %change in bottom distance
            netbdt = NaN(maxNumEvnt,1); %rate of change in bottom distance
            deltaT = NaN(maxNumEvnt,1); %event length (timesteps*dt)
            bdmed = NaN(maxNumEvnt,1); %median bed level
            umag = NaN(maxNumEvnt,1); %velocity magnitude
            usqd = NaN(maxNumEvnt,1); %velocity squared
            umed = NaN(maxNumEvnt,1); %median velocity
            phase = NaN(maxNumEvnt,1); %phase (deg) from wavelet
            sigh = NaN(maxNumEvnt,1); %significant wave height from aquadopp
            depth = NaN(maxNumEvnt,1); %water depth from aquadopp
            orbwv = NaN(maxNumEvnt,1); %orbital wave vel from aquadopp
            tbd = NaN(maxNumEvnt,1); %bed shear stress
            epd = NaN(maxNumEvnt,1); %turbulent dissipation
            disp('Determining the max length of each event...')
            t3 = tic;
            for k = 1:maxNumEvnt
                %find each event, then find the longest axis of that event. This number
                %is going to be the event length
                [rw,co] = find(eventGroups == k); %find indices of a single event
                unrw = unique(rw);
                eventLength = zeros(length(unrw),1);
                for kk = 1:length(unrw)           %figure out which row is longest
                    eventLength(kk) = numel(co(rw == unrw(kk)));
                end
                [maxEvent,maxID] = max(eventLength);
                gid = co(rw == unrw(maxID));
                
                %Continue script as in EventDuration_v2.mat
                nevent(k) = numel(gid); %# timesteps in event
                deltabd(k) = bds(gid(end))-bds(gid(1)); %net bed change
                deltaT(k) = maxEvent*0.1; %time in seconds
                netbdt(k) = deltabd(k)/deltaT(k); %net rate of bed change
                bdmed(k) = nanmedian(bds(gid)); %median bed level during event
                umag(k) = nanmean(abs(x10(gid))); %average velocity magnitude over event
                umed(k) = nanmedian(x10(gid)); %median velocity
                usqd(k) = nanmean(x10(gid).^2); %velocity squared over event
                phase(k) = angle(nanmean(wvlt.(wtdir{j}).Wxy(unrw(maxID),gid)))*(180/pi); %phase in deg.
                sigh(k) = nanmean(Hs10(gid)); %sig wave height
                depth(k) = nanmean(H10(gid)); %water depth
                orbwv(k) = nanmean(Om10(gid)); %orbital wave velocity
                tbd(k) = nanmean(tb10(gid)); %bed shear stress
                epd(k) = nanmean(ep10(gid)); %turbulence
            end
            fprintf('Event length extraction completed in %0.2f minutes\n',toc(t3)/60)            
            %% Filter by wave/ig
            wvband = find(deltaT<30); %less than 30 seconds
            igband = find(deltaT>=30); %greater than 30 seconds
            wvevent = deltaT(wvband,:);
            igevent = deltaT(igband,:);
            wvdbd = deltabd(wvband,:);
            igdbd = deltabd(igband,:);
            wvnbd = netbdt(wvband,:);
            ignbd = netbdt(igband,:);
            wvbdm = bdmed(wvband,:);
            igbdm = bdmed(igband,:);
            wvumd = umed(wvband,:);
            igumd = umed(igband,:);
            wvusq = usqd(wvband,:);
            igusq = usqd(igband,:);
            wvumg = umag(wvband,:);
            igumg = umag(igband,:);
            wvphs = phase(wvband,:);
            igphs = phase(igband,:);
            wvsgh = sigh(wvband,:);
            igsgh = sigh(igband,:);
            wvdep = depth(wvband,:);
            igdep = depth(igband,:);
            wvorb = orbwv(wvband,:);
            igorb = orbwv(igband,:);
            wvtbd = tbd(wvband,:);
            igtbd = tbd(igband,:);
            wvepd = epd(wvband,:);
            igepd = epd(igband,:);
            %% Save out to structure
            eventd.(dfn{rvp(vpid(j))}).wave.bdst = bds(find(bds~=0,1,'first')); %first bottom elevation estimate
            eventd.(dfn{rvp(vpid(j))}).wave.eventl = wvevent; %event length
            eventd.(dfn{rvp(vpid(j))}).wave.deltbd = wvdbd; %net change in bed level in event
            eventd.(dfn{rvp(vpid(j))}).wave.netbdt = wvnbd; %rate of net change in bed level in event
            eventd.(dfn{rvp(vpid(j))}).wave.bdmed = wvbdm; %median bed level change in event
            eventd.(dfn{rvp(vpid(j))}).wave.umed = wvumd; %median velocity in event
            eventd.(dfn{rvp(vpid(j))}).wave.usqd = wvusq; %velocity squared in event
            eventd.(dfn{rvp(vpid(j))}).wave.umag = wvumg; %velocity magnitude in x-direction
            eventd.(dfn{rvp(vpid(j))}).wave.phase = wvphs; %avg. phase of wavelet coherence
            eventd.(dfn{rvp(vpid(j))}).wave.sigh = wvsgh; %significant wave height in event
            eventd.(dfn{rvp(vpid(j))}).wave.depth = wvdep; %water depth in event
            eventd.(dfn{rvp(vpid(j))}).wave.orbwv = wvorb; %wave orbital velocity in event
            eventd.(dfn{rvp(vpid(j))}).wave.taub = wvtbd; %bed shear stress in event
            eventd.(dfn{rvp(vpid(j))}).wave.eps = wvepd; %turbulence in event
            eventd.(dfn{rvp(vpid(j))}).ig.bdst = bds(find(bds~=0,1,'first'));
            eventd.(dfn{rvp(vpid(j))}).ig.eventl = igevent;
            eventd.(dfn{rvp(vpid(j))}).ig.deltbd = igdbd;
            eventd.(dfn{rvp(vpid(j))}).ig.netbdt = ignbd;
            eventd.(dfn{rvp(vpid(j))}).ig.bdmed = igbdm;
            eventd.(dfn{rvp(vpid(j))}).ig.umed = igumd;
            eventd.(dfn{rvp(vpid(j))}).ig.usqd = igusq;
            eventd.(dfn{rvp(vpid(j))}).ig.umag = igumg;
            eventd.(dfn{rvp(vpid(j))}).ig.phase = igphs;
            eventd.(dfn{rvp(vpid(j))}).ig.sigh = igsgh;
            eventd.(dfn{rvp(vpid(j))}).ig.depth = igdep;
            eventd.(dfn{rvp(vpid(j))}).ig.orbwv = igorb;
            eventd.(dfn{rvp(vpid(j))}).ig.taub = igtbd;
            eventd.(dfn{rvp(vpid(j))}).ig.eps = igepd;
            clear dat wvlt tke tau
        end
        fprintf('File processing completed in %0.2f minutes\n',toc(t1)/60)
        save([npath 'Wavelet\' folders{ii} '\' 'waveIGevents'],'eventd','-v7.3')
    end
end