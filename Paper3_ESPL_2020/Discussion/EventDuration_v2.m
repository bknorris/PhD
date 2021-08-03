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
rdir = 'g:\GradSchool\DataAnalysis\Paper3\ExperimentalDesign\';
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
	npath = ['f:\' ypath{i} dpath];
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
        eventd = struct();
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
            %% Filter by movement events
            %find high coherence events (sig95 > 0.9) per bandwidth
            threshold = 0.9;
            [m,n] = size(wvlt.(wtdir{j}).sig95);
            t = wvlt.(wtdir{j}).t;
            events = zeros(m,n);
            for k = 1:m
                events(k,:) = wvlt.(wtdir{j}).sig95(k,:) >= threshold;
            end
            %filter by coi
            zid = zeros(m,n);
            for k = 1:n
                zid(:,k) = wvlt.(wtdir{j}).period <= wvlt.x.coi(k);
            end
            events = events.*zid;
            events(events==0) = NaN;
            %count length of gap using cmgidgaps
            nevent = NaN(m,n);lgap = zeros(m,1);
            bdvn = NaN(m,n); %normalized bottom variance
            deltabd = NaN(m,n); %change in bottom distance
            bdmed = NaN(m,n); %median bed level
            umag = NaN(m,n); %velocity magnitude
            usqd = NaN(m,n); %velocity squared
            umed = NaN(m,n); %median velocity
            phase = NaN(m,n); %phase (deg) from wavelet
            sigh = NaN(m,n); %significant wave height from aquadopp
            depth = NaN(m,n); %water depth from aquadopp
            orbwv = NaN(m,n); %orbital wave vel from aquadopp
            tbd = NaN(m,n); %bed shear stress
            epd = NaN(m,n); %turbulent dissipation
            for k = 1:m
                [gap,last,first,~] = cmgidgaps(events(k,:));
                for kk = 1:gap-1
                    %calc sigma_bd/deltaTU for sig95 > 0.9 events
                    gid = first(kk)+1:last(kk+1)-1;
                    U = nanmean(abs(x10(gid))); %average velocity magnitude over event
                    umag(k,kk) = U;
                    usqd(k,kk) = nanmean(x10(gid).^2); %velocity squared over event
                    deltaT = (t(gid(end))-t(gid(1)))*60; %time in seconds
                    sigbd = var(bds(gid)); %variance of bed over event
                    nevent(k,kk) = numel(gid); %# timesteps in event
                    if isnan(sigbd/(U*deltaT))
                        bdvn(k,kk) = 0;
                    else
                        bdvn(k,kk) = sigbd/(U*deltaT); %normalized bottom variance
                    end
                    deltabd(k,kk) = bds(gid(end))-bds(gid(1)); %net bed change
                    bdmed(k,kk) = nanmedian(bds(gid)); %median bed level
                    umed(k,kk) = nanmedian(x10(gid)); %median velocity
                    phase(k,kk) = angle(nanmean(wvlt.(wtdir{j}).Wxy(k,gid)))*(180/pi); %phase in deg.
                    sigh(k,kk) = nanmean(Hs10(gid)); %sig wave height
                    depth(k,kk) = nanmean(H10(gid)); %water depth
                    orbwv(k,kk) = nanmean(Om10(gid)); %orbital wave velocity
                    tbd(k,kk) = nanmean(tb10(gid)); %bed shear stress
                    epd(k,kk) = nanmean(ep10(gid));
                end
                lgap(k) = gap;
            end
            %get rid of extra NaNs
            nevent = nevent(:,1:max(lgap)-1); 
            bdvn = bdvn(:,1:max(lgap)-1);
            deltabd = deltabd(:,1:max(lgap)-1);
            bdmed = bdmed(:,1:max(lgap)-1);
            umag = umag(:,1:max(lgap)-1);
            usqd = usqd(:,1:max(lgap)-1);
            umed = umed(:,1:max(lgap)-1);
            phase = phase(:,1:max(lgap)-1);
            sigh = sigh(:,1:max(lgap)-1);
            depth = depth(:,1:max(lgap)-1);
            orbwv = orbwv(:,1:max(lgap)-1);
            tbd = tbd(:,1:max(lgap)-1);
            epd = epd(:,1:max(lgap)-1);
            dt = t(2)-t(1);
            ntime = dt.*nevent;
            %% Filter by wave-band vs. IG
            period = wvlt.(wtdir{j}).period;
            wvband = find(period<=0.125); 
            igband = find(period>1);
            wvevent = ntime(wvband,:);
            igevent = ntime(igband,:);
            wvbdvn = bdvn(wvband,:);
            igbdvn = bdvn(igband,:);
            wvdbd = deltabd(wvband,:);
            igdbd = deltabd(igband,:);
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
            eventd.(dfn{rvp(vpid(j))}).wave.bdvar = wvbdvn; %normalized bottom variance
            eventd.(dfn{rvp(vpid(j))}).wave.deltbd = wvdbd; %net change in bed level in event
            eventd.(dfn{rvp(vpid(j))}).wave.bdmed = wvbdm; %median bed level change in event
            eventd.(dfn{rvp(vpid(j))}).wave.umed = wvumd; %median velocity in event
            eventd.(dfn{rvp(vpid(j))}).wave.usqd = wvusq; %velocity squared in event
            eventd.(dfn{rvp(vpid(j))}).wave.umag = wvumg; %velocity magnitude in x-direction
            eventd.(dfn{rvp(vpid(j))}).wave.phase = wvphs; %avg. phase of wavelet coherence
            eventd.(dfn{rvp(vpid(j))}).wave.period = period(wvband); %period in waveband
            eventd.(dfn{rvp(vpid(j))}).wave.sigh = wvsgh; %significant wave height in event
            eventd.(dfn{rvp(vpid(j))}).wave.depth = wvdep; %water depth in event
            eventd.(dfn{rvp(vpid(j))}).wave.orbwv = wvorb; %wave orbital velocity in event
            eventd.(dfn{rvp(vpid(j))}).wave.taub = wvtbd; %bed shear stress in event
            eventd.(dfn{rvp(vpid(j))}).wave.eps = wvepd; %turbulence in event
            eventd.(dfn{rvp(vpid(j))}).ig.bdst = bds(find(bds~=0,1,'first'));
            eventd.(dfn{rvp(vpid(j))}).ig.eventl = igevent;
            eventd.(dfn{rvp(vpid(j))}).ig.bdvar = igbdvn;
            eventd.(dfn{rvp(vpid(j))}).ig.deltbd = igdbd;
            eventd.(dfn{rvp(vpid(j))}).ig.bdmed = igbdm;
            eventd.(dfn{rvp(vpid(j))}).ig.umed = igumd;
            eventd.(dfn{rvp(vpid(j))}).ig.usqd = igusq;
            eventd.(dfn{rvp(vpid(j))}).ig.umag = igumg;
            eventd.(dfn{rvp(vpid(j))}).ig.phase = igphs;
            eventd.(dfn{rvp(vpid(j))}).ig.period = period(igband);
            eventd.(dfn{rvp(vpid(j))}).ig.sigh = igsgh;
            eventd.(dfn{rvp(vpid(j))}).ig.depth = igdep;
            eventd.(dfn{rvp(vpid(j))}).ig.orbwv = igorb;
            eventd.(dfn{rvp(vpid(j))}).ig.taub = igtbd;
            eventd.(dfn{rvp(vpid(j))}).ig.eps = igepd;
            clear dat wvlt tke tau
        end
        save([npath 'Wavelet\' folders{ii} '\' 'waveIGevents'],'eventd','-v7.3')
    end
end