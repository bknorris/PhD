%need to fix how an event is defined in the wavelet plots.
clear
load('F:\Mekong_W2015\DataAnalysis\Paper3\Wavelet\08-03-15\vpro2_wvlt.mat')
%% Filter by movement events
%find high coherence events (sig95 > 0.9) per bandwidth
threshold = 0.9;
[m,n] = size(wvlt.x.sig95);
t = wvlt.x.t;
events = zeros(m,n);
for k = 1:m
    events(k,:) = wvlt.x.sig95(k,:) >= threshold;
end
%filter by coi
zid = zeros(m,n);
for k = 1:n
    zid(:,k) = wvlt.x.period <= wvlt.x.coi(k);
end
events = events.*zid;
%% do some 'image' processing; WARNING: Requires image processing toolbox!!!
eventGroups = bwlabel(events,8);
maxNumEvnt = max(max(eventGroups));
%initialize variables to save out
nevent = NaN(maxNumEvnt,1); %# timesteps in event
deltabd = NaN(maxNumEvnt,1); %change in bottom distance
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
t1 = tic;
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
%     deltabd(k) = bds(gid(end))-bds(gid(1)); %net bed change
    deltaT(k) = maxEvent*dt; %time in seconds
%     bdmed(k) = nanmedian(bds(gid)); %median bed level during event
%     umag(k) = nanmean(abs(x10(gid))); %average velocity magnitude over event
%     umed(k) = nanmedian(x10(gid)); %median velocity
%     usqd(k) = nanmean(x10(gid).^2); %velocity squared over event
%     phase(k) = angle(nanmean(wvlt.(wtdir{j}).Wxy(maxID,gid)))*(180/pi); %phase in deg.
%     sigh(k) = nanmean(Hs10(gid)); %sig wave height
%     depth(k) = nanmean(H10(gid)); %water depth
%     orbwv(k) = nanmean(Om10(gid)); %orbital wave velocity
%     tbd(k) = nanmean(tb10(gid)); %bed shear stress
%     epd(k) = nanmean(ep10(gid)); %turbulence
end
fprintf('Event length extraction completed in %0.2f minutes\n',toc(t1)/60)