%Check all wavelet phases for comment-response, ESPL
%Norris, 2020
clear, close all

%define working paths
dpath = '\DataAnalysis\Paper3\';
ypath = 'Mekong_W2015';
%load run file to tell program which files & VPs to load
rdir = 'd:\MemoryStick\GradSchool\DataAnalysis\Paper3\ExperimentalDesign\';
fid = fopen([rdir 'AutoMLRrunfile_flood.csv']);
rfile = textscan(fid,'%s%s%s%n%n%n%s%s','delimiter',',');
rdate = rfile{1};
dfn = {'vpro1';'vpro2';'vpro3'};
npath = ['d:\' ypath dpath];
folders = dir([npath '\Wavelet\']);
folders = {folders(3:end).name};
for ii = 1:length(folders)
    if any(strcmp(folders{ii},rdate))
        vpid = strcmp(folders{ii},rdate);vpid = find(vpid);
        disp(['Loading files from ' npath 'Wavelet\' folders{ii} '\'])
        wfile  = dir([npath 'Wavelet\' folders{ii} '\','*.mat']);
    else
        continue
    end
    for j = 2:length(vpid)
        %% Load files
        disp(['Processing ' dfn{j}])
        load([npath 'Wavelet\' folders{ii} '\' wfile(j).name]);
        sig95 = wvlt.x.sig95;
        Rsq = wvlt.x.Rsq;
        Wxy = wvlt.x.Wxy;
        coi = wvlt.x.coi;
        t = wvlt.x.t;
        period = wvlt.x.period;
        dt = wvlt.x.dt;
        
        %find high coherence events (sig95 > 0.9) per bandwidth
        threshold = 0.9;
        [m,n] = size(sig95);
        events = zeros(m,n);
        for k = 1:m
            events(k,:) = sig95(k,:) >= threshold;
        end
        %filter by coi
        zid = zeros(m,n);
        for k = 1:n
            zid(:,k) = wvlt.x.period <= coi(k);
        end
        events = events.*zid;
        %% do some 'image' processing; WARNING: Requires image processing toolbox!!!
        eventGroups = bwlabel(events,8);
        maxNumEvnt = max(max(eventGroups));
        events = eventGroups~=0; %just the events
        aWxy = angle(Wxy.*events); %phase in radians
        avgaWxy = NaN(m,n);
        for k = 1:maxNumEvnt
            %find each event, then average angle in each event
            aidx = find(eventGroups == k); %find indices of a single event
            avgaWxy(aidx) = mean(aWxy(aidx));
            %           avgaWxy(aidx) = aWxy(aidx);
        end
        %% Plot results on a hist
        figure
        igband = find(period>1 & period < 35);
        aevents = avgaWxy(igband,:);aevents = aevents(~isnan(aevents)); %find phases in events within the IG band
        hist(rad2deg(aevents),-180:45:180)
        set(gca,'xlim',[-200 200])
        xlabel('Phase angle (deg)')
        ylabel('Counts')
        titletext = ['Phase analysis - ' folders{ii} ' ' dfn{j}];
        title(titletext)
    end
    pause()
    close all
end
        
