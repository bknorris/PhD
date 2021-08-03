%Extract phase from the wavelet plots for 03/09/15 (flood-ebb)
%
clear, close all
%
%% First plot the bed level time series from 03/09
ddir = 'd:\Mekong_W2015\DataAnalysis\Paper3\BottomTrack\';
load([ddir 'BDtrace_floodebb.mat'])
f1 = figure(1);
set(f1,'PaperOrientation','portrait',...
    'position',[400 100   1000  350],...
    'renderer','painters');
names = {'vp1';'vp2';'vp3'};
% symb = {'o';'d';'p'};
lines = {'-';'-';'-'};
tstep = datenum(0,0,0,0,60,0);
hight = datenum(2015,03,09,04,33,00);
cl = flipud([0.1 0.1 0.1;0.5 0.5 0.5;0.7 0.7 0.7]);
pl = zeros(3,1);
space = 300;
for i = 1:3
    time = data.(names{i}).time;
    lvl = data.(names{i}).ble.*1000;
    if i == 1    
        plot(time,zeros(length(time),1),'--k',...
            'linewidth',1.5),hold on
        plot(hight*ones(10,1),linspace(-100,100,10),'-k',...
            'linewidth',1.5)
    end
    pl(i) = plot(time,lvl,lines{i},...
        'color',cl(i,:),'linewidth',1.5);
    lc = 1:space:length(time);
end
leg = legend(pl,{'Seaward';'Middle';'Landward'});
datetickzoom('x','dd HH:MM:SS','keepticks','keeplimits')
xlabel('Time on 09-03-15 [HH:MM]')
ylabel('Bed Level [mm]')

%% Run the wavelet analysis
%Define working paths
dpath = 'd:\Mekong_W2015\DataAnalysis\Paper3\Wavelet\09-03-15\';
figdir = 'd:\Mekong_W2015\Figures\Paper3\Wavelet\';
fname = dir([dpath '*_wvlt.mat']);fname = {fname.name};
%We need VPRO1 and VPRO2 for this analysis
insts = {'vpro1';'vpro2'};
for i = 1:2
    fid = find(contains(fname,insts{i}));
    disp(['Loading ' dpath fname{fid}])
    load([dpath fname{fid}])
    %% Filter by movement events
    %crop all variables to flood-high and high-ebb for i = 1 and i = 2,
    %respectively
    ht = floor(length(wvlt.x.t)/2);
    if i == 1 %this is the flood tide
        tidx = 1:ht;
    else %this is the ebb tide
        tidx = ht:length(wvlt.x.t);
    end
    sig95 = wvlt.x.sig95(:,tidx);
    Rsq = wvlt.x.Rsq(:,tidx);
    Wxy = wvlt.x.Wxy(:,tidx);
    coi = wvlt.x.coi(tidx);
    t = wvlt.x.t(tidx);
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
            
    %Plot figure as a reference
    cb = brewermap(100,'*RdYlBu');
    f2 = figure(2);
    plotwtc(Rsq,period,coi,sig95,t,Wxy,dt)
    colormap(gca,cb)
    xlabel('Time (minutes)'),ylabel('Period (minutes)')
    cbl = colorbar('peer',gca,'Tag','colorbar1');
    ylabel(cbl,'Magnitude Squared Coherence')
    title(['Wavelet, 09/03/15 - ' fname{fid}(1:5)])
    
    igband = find(period>1 & period < 35);
    aevents = avgaWxy(igband,:);aevents = aevents(~isnan(aevents)); %find phases in events within the IG band
    figure(3);
    imagesc(t,period,avgaWxy)
    figure(4)
    imagesc(t,period(period>1&period<35),rad2deg(avgaWxy(period>1&period<35,:)))
    f3 = figure(5);
    hist(rad2deg(aevents),-180:45:180)
    set(gca,'xlim',[-200 200])
    xlabel('Phase angle (deg)')
    ylabel('Counts')
    title(['Phase Angle from Coherent "events", 09/03/15 - ' fname{i}(1:5)])
    
    bins = -180:45:180;
    h = hist(rad2deg(aevents),bins);
    pctpos = (sum(h(1:4))/sum(h))*100
    %Save figures
%     prettyfigures('text',10,'labels',12,'box',1)
%     if i == 1
%         export_fig(f1,[figdir 'Bedlevels_090315'],'-png')
%     end
%     export_fig(f2,[figdir fname{fid}(1:5) '_Wavelet'],'-png')
%     export_fig(f3,[figdir fname{fid}(1:5) '_phaseHist'],'-png')
end

