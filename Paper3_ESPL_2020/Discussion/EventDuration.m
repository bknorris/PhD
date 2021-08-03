%Try to extract the event duration from individual frequency bands of the
%wavelet analysis (i.e., how long a high-coherence event lasts). Plot this
%against the normalized bed level elevation variance (e.g., Staudt et al.
%2017)
clear
%define working paths
dpath = '\DataAnalysis\Paper3\';
ypath = {'Mekong_F2014';'Mekong_W2015'};
%load run file to tell program which files & VPs to load
rdir = 'e:\GradSchool\DataAnalysis\Paper3\ExperimentalDesign\';
fid = fopen([rdir 'AutoEventDrunfile.csv']);
rfile = textscan(fid,'%s%s%s%n%n%n%s','delimiter',',');
rdate = rfile{1};rstart = rfile{2};%dates corresponding to folder names; start time for crop
rstop = rfile{3};rvp = rfile{4};  %stop time for crop; vp #
rbin5 = rfile{5};rhead = rfile{6}; %use bins 1-5 for averaging; VP heading for cross & along-shore rot
wtdir = rfile{7};                  %wavelet direction ('x' or 'y')
start = datenum(strcat(rdate,{' '},rstart),'dd-mm-yy HH:MM:SS');
stop = datenum(strcat(rdate,{' '},rstop),'dd-mm-yy HH:MM:SS');
dfn = {'vpro1';'vpro2';'vpro3'};
for i = 1:2
	npath = ['g:\' ypath{i} dpath];
    folders = dir([npath '\VPs\']);
    folders = {folders(3:end).name};
    for ii = 1:length(folders)
        if any(strcmp(folders{ii},rdate))
            vpid = strcmp(folders{ii},rdate);vpid = find(vpid);
            disp(['Loading files from ' npath 'VPs\' folders{ii} '\'])
            disp(['Loading files from ' npath 'BottomTrack\' folders{ii} '\'])
            vfile = dir([npath 'VPs\' folders{ii} '\','*_Vels.mat']);
            wfile  = dir([npath 'Wavelet\' folders{ii} '\','*.mat']); 
            bfile = dir([npath 'BottomTrack\' folders{ii} '\','*_bdtrace.mat']);
            bd = load([npath 'BottomTrack\' folders{ii} '\' bfile.name]);
        else
            continue
        end
        %just load the variables you need, saves memory
        eventd = struct();
        for j = 1:length(vpid)
            disp(['Processing ' dfn{rvp(vpid(j))}])
            dat = load([npath 'VPs\' folders{ii} '\' vfile.name],dfn{rvp(vpid(j))});
            load([npath 'Wavelet\' folders{ii} '\' wfile(j).name]);
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
            umed = NaN(m,n); %median velocity
            phase = NaN(m,n); %phase (deg) from wavelet
            for k = 1:m
                [gap,last,first,~] = cmgidgaps(events(k,:));
                for kk = 1:gap-1
                    %calc sigma_bd/deltaTU for sig95 > 0.9 events
                    gid = first(kk)+1:last(kk+1)-1;
                    U = nanmean(abs(x10(gid))); %average velocity magnitude over event
                    umag(k,kk) = U;
                    deltaT = (t(gid(end))-t(gid(1)))*60; %time in seconds
                    sigbd = var(bds(gid)); %variance of bed over event
                    nevent(k,kk) = numel(gid); %# timesteps in event
                    bdvn(k,kk) = sigbd/(U*deltaT); %normalized bottom variance
                    deltabd(k,kk) = bds(gid(end))-bds(gid(1)); %net bed change
                    bdmed(k,kk) = nanmedian(bds(gid)); %median bed level
                    umed(k,kk) = nanmedian(x10(gid)); %median velocity
                    phase(k,kk) = angle(nanmean(wvlt.(wtdir{j}).Wxy(k,gid)))*(180/pi); %phase in deg.
                end
                lgap(k) = gap;
            end
            %get rid of extra NaNs
            nevent = nevent(:,1:max(lgap)-1); 
            bdvn = bdvn(:,1:max(lgap)-1);
            deltabd = deltabd(:,1:max(lgap)-1);
            bdmed = bdmed(:,1:max(lgap)-1);
            umag = umag(:,1:max(lgap)-1);
            umed = umed(:,1:max(lgap)-1);
            phase = phase(:,1:max(lgap)-1);
            dt = t(2)-t(1);
            ntime = dt.*nevent;
            %now filter by wave-band vs. IG
            period = wvlt.(wtdir{j}).period;
            wvband = find(period<=0.125); %7.5 second period (likely an overestimate)
            igband = find(period>0.125);
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
            wvumg = umag(wvband,:);
            igumg = umag(igband,:);
            wvphs = phase(wvband,:);
            igphs = phase(igband,:);
            %save out to structure
            eventd.(dfn{rvp(vpid(j))}).wave.bdst = bds(find(bds~=0,1,'first')); %first bottom elevation estimate
            eventd.(dfn{rvp(vpid(j))}).wave.eventl = wvevent; %event length
            eventd.(dfn{rvp(vpid(j))}).wave.bdvar = wvbdvn; %normalized bottom variance
            eventd.(dfn{rvp(vpid(j))}).wave.deltbd = wvdbd; %net change in bed level in event
            eventd.(dfn{rvp(vpid(j))}).wave.bdmed = wvbdm; %median bed level change in event
            eventd.(dfn{rvp(vpid(j))}).wave.umed = wvumd; %median velocity in event
            eventd.(dfn{rvp(vpid(j))}).wave.umag = wvumg; %velocity magnitude in x-direction
            eventd.(dfn{rvp(vpid(j))}).wave.phase = wvphs; %avg. phase of wavelet coherence
            eventd.(dfn{rvp(vpid(j))}).wave.period = period(wvband); %period in waveband
            eventd.(dfn{rvp(vpid(j))}).ig.bdst = bds(find(bds~=0,1,'first'));
            eventd.(dfn{rvp(vpid(j))}).ig.eventl = igevent;
            eventd.(dfn{rvp(vpid(j))}).ig.bdvar = igbdvn;
            eventd.(dfn{rvp(vpid(j))}).ig.deltbd = igdbd;
            eventd.(dfn{rvp(vpid(j))}).ig.bdmed = igbdm;
            eventd.(dfn{rvp(vpid(j))}).ig.umed = igumd;
            eventd.(dfn{rvp(vpid(j))}).ig.umag = igumg;
            eventd.(dfn{rvp(vpid(j))}).ig.phase = igphs;
            eventd.(dfn{rvp(vpid(j))}).ig.period = period(igband);
            clear dat wvlt
        end
        save([npath 'Wavelet\' folders{ii} '\' 'waveIGevents'],'eventd','-v7.3')
    end
end