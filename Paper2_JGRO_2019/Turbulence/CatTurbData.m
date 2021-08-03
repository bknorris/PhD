%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%We're going to plot TKE dissipation rates (averaged over 10 minute
%windows) against the vegetation statistics (phi) to see if there is a
%trend in TKE with increasing vegetation density. In this case, plot PHI by
%CROSS-SHORE DISTANCE (X) and compare this to TKE dissipation rates by
%distance.

%This script is to process the velocity data into a sensible structure,
%which will also include wave data (akin to the VegDensityTKE.m script).
%Distance values are located in the dat structure from CatVegGeometry.m
%The actual plots will be made by another script:
%
%
%Velocity data records(et al) will be divided into three time periods: low,
%mid and high tide as evenly as possible. THese divisions can be recombined
%for the eventual TKE/Phi/Distance figure.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%% Load the 2014 Data %%%%
veldir = 'd:\Projects\Mekong_F2014\DataAnalysis\Paper2\AveragedVelocities\';
wvdir = 'd:\Projects\Mekong_F2014\DataAnalysis\Paper2\Spectra\';
tkedir = 'd:\Projects\Mekong_F2014\DataAnalysis\Paper2\TKE\';
velfiles = dir(veldir);velfiles = {velfiles.name};
tkefiles = dir(tkedir);tkefiles = {tkefiles.name};
wvfiles = dir(wvdir);wvfiles = {wvfiles.name};
%%%%extract TKE estimates from TKE files%%%%

horiz = [1 1 0 0 1 0 1 0];
ebb = [0 1 0 1 0 0 0 0];                                                    %experiments where the ebb tide was recorded
eps = zeros(24,3);
avgs = zeros(24,3);
stdev = zeros(24,3);
wvrms = zeros(24,1);
depth = zeros(24,1);
name = cell(8,1);
e = zeros(3,3);avg = zeros(3,3);wrms = zeros(3,1);
tide = zeros(3,1);stdv = zeros(3,3);
c = 1;cc = 1;                                                               %index counters
for i = 3:length(tkefiles)
    load([tkedir tkefiles{i}])
    load([veldir velfiles{i}])
    load([wvdir wvfiles{i}])
    name{c} = regexprep(velfiles{i},'Vave.mat','');
    
    %Use times to determine the number of samples to average over by
    %dividing the timeseries in thirds
    start = Stat.vpro1.time(1);stop = Stat.vpro1.time(end);
    [~,~,~,hr,mi,~] = datevec(stop-start);
    time = (hr*60) + mi;                                                    %calculate minutes elapsed between start/stop time of experiment
    time = datenum(0,0,0,0,floor(time/10),0);
    if ebb(c) == 0
        timee = start:time:stop;timee(end) = stop;                          %if the data is from an ebb tide, flip the time record around so the first measurement is lower tidal elevation
    elseif ebb(c) == 1
        timee = fliplr(start:time:stop);timee(1) = stop;
    end
    for j = 1:length(timee)-1
        if ebb(c) == 0
            id = find(Stat.vpro1.time >= timee(j) & Stat.vpro1.time <= timee(j+1));
            id2 = find(wvstats.time >= timee(j) & wvstats.time <= timee(j+1));
        elseif ebb(c) == 1
            id = find(Stat.vpro1.time <= timee(j) & Stat.vpro1.time >= timee(j+1));
            id2 = find(wvstats.time <= timee(j) & wvstats.time >= timee(j+1));
        end
        if horiz(c) == 1                                                    %if the instrument is horizontal, use bin 15 for TKE estimates
            tk(1) = nanmean((Stat.vpro1.beam1.E(15,id)+Stat.vpro1.beam3.E(15,id))./2);
            av(1) = nanmean(Avgs.vpro1.x(15,id));
            sd(1) = std((Stat.vpro1.beam1.E(15,id)+Stat.vpro1.beam3.E(15,id))./2);
            tk(2) = nanmean((Stat.vpro2.beam1.E(15,id)+Stat.vpro2.beam3.E(15,id))./2);
            av(2) = nanmean(Avgs.vpro2.x(15,id));
            sd(2) = std((Stat.vpro2.beam1.E(15,id)+Stat.vpro2.beam3.E(15,id))./2);
            tk(3) = nanmean((Stat.vpro3.beam1.E(15,id)+Stat.vpro3.beam3.E(15,id))./2);
            av(3) = nanmean(Avgs.vpro3.x(15,id));
            sd(3) = std((Stat.vpro3.beam1.E(15,id)+Stat.vpro3.beam3.E(15,id))./2);
            hrms = nanmean(wvstats.hrmsp(id2));
            dep = nanmean(wvstats.depth(id2));
        elseif horiz(c) == 0
            tk(1) = nanmean((Stat.vpro1.beam1.E(5,id)+Stat.vpro1.beam3.E(5,id))./2);
            av(1) = nanmean(Avgs.vpro1.x(5,id));
            sd(1) = std((Stat.vpro1.beam1.E(5,id)+Stat.vpro1.beam3.E(5,id))./2);
            tk(2) = nanmean((Stat.vpro2.beam1.E(5,id)+Stat.vpro2.beam3.E(5,id))./2);
            av(2) = nanmean(Avgs.vpro2.x(5,id));
            sd(2) = std((Stat.vpro2.beam1.E(5,id)+Stat.vpro2.beam3.E(5,id))./2);
            tk(3) = nanmean((Stat.vpro3.beam1.E(5,id)+Stat.vpro3.beam3.E(5,id))./2);
            av(3) = nanmean(Avgs.vpro3.x(5,id));
            sd(3) = std((Stat.vpro3.beam1.E(5,id)+Stat.vpro3.beam3.E(5,id))./2);
            hrms = nanmean(wvstats.hrmsp(id2));
            dep = nanmean(wvstats.depth(id2));
        end
        e(j,:) = tk;
        avg(j,:) = av;
        wrms(j,:) = hrms;
        tide(j,:) = dep;
        stdv(j,:) = sd;
    end
    eps(cc:cc+2,:) = e;
    avgs(cc:cc+2,:) = abs(avg);
    wvrms(cc:cc+2,:) = wrms;
    depth(cc:cc+2,:) = tide;
    stdev(cc:cc+2,:) = stdv;
    c = c+1;
    cc = cc+3;
end
%rearrange the data in the order of deployment (Q2-Q7).
E = NaN(24,3);
E(1:3,:) = eps(22:24,:);E(4:6,:) = eps(16:18,:);E(7:9,:) = eps(19:21,:);
E(10:12,:) = eps(1:3,:);E(13:15,:) = eps(13:15,:);E(16:18,:) = eps(7:9,:);
E(19:21,:) = eps(10:12,:);E(22:24,:) = eps(4:6,:);
Avg = NaN(24,3);
Avg(1:3,:) = avgs(22:24,:);Avg(4:6,:) = avgs(16:18,:);Avg(7:9,:) = avgs(19:21,:);
Avg(10:12,:) = avgs(1:3,:);Avg(13:15,:) = avgs(13:15,:);Avg(16:18,:) = avgs(7:9,:);
Avg(19:21,:) = avgs(10:12,:);Avg(22:24,:) = avgs(4:6,:);
Depth = NaN(24,1);
Depth(1:3,:) = depth(22:24,:);Depth(4,:) = depth(17,:);Depth(5:6,:) = NaN;
Depth(7:9,:) = depth(19:21,:);Depth(10:12,:) = depth(1:3,:);
Depth(13:15,:) = depth(13:15,:);Depth(16:18,:) = depth(7:9,:);
Depth(19:21,:) = depth(10:12,:);Depth(22:24,:) = depth(4:6,:);
Depth = repmat(Depth,1,3);
Wrms = NaN(24,1);
Wrms(1:3,:) = wvrms(22:24,:);Wrms(4,:) = wvrms(17,:);Wrms(5:6,:) = NaN;
Wrms(7:9,:) = wvrms(19:21,:);Wrms(10:12,:) = wvrms(1:3,:);
Wrms(13:15,:) = wvrms(13:15,:);Wrms(16:18,:) = wvrms(7:9,:);
Wrms(19:21,:) = wvrms(10:12,:);Wrms(22:24,:) = wvrms(4:6,:);
Wrms = repmat(Wrms,1,3);
Stdev = zeros(24,3);
Stdev(1:3,:) = stdev(22:24,:);Stdev(4,:) = stdev(17,:);Stdev(5:6,:) = NaN;
Stdev(7:9,:) = stdev(19:21,:);Stdev(10:12,:) = stdev(1:3,:);
Stdev(13:15,:) = stdev(13:15,:);Stdev(16:18,:) = stdev(7:9,:);
Stdev(19:21,:) = stdev(10:12,:);Stdev(22:24,:) = stdev(4:6,:);
quad = cell(24,3);
quad(1:3,:) = name(8);quad(4:6,:) = name(6);
quad(7:9,:) = name(7);quad(10:12,:) = name(1);
quad(13:15,:) = name(5);quad(16:18,:) = name(3);
quad(19:21,:) = name(4);quad(22:24,:) = name(2);
veldat.four.quad = quad;veldat.four.E = E;veldat.four.Estd = Stdev;
veldat.four.avg = Avg;veldat.four.depth = Depth;veldat.four.wrms = Wrms;

clearvars -except veldat vegdat

%%%% Load the 2015 Data %%%%
veldir = 'd:\Projects\Mekong_W2015\DataAnalysis\Paper2\AveragedVelocities\';
tkedir = 'd:\Projects\Mekong_W2015\DataAnalysis\Paper2\TKE\';
wvdir = 'd:\Projects\Mekong_W2015\DataAnalysis\Paper2\Spectra\';
velfiles = dir(veldir);velfiles = {velfiles.name};
tkefiles = dir(tkedir);tkefiles = {tkefiles.name};
wvfiles = dir(wvdir);wvfiles = {wvfiles.name};

eps = zeros(27,3);
avgs = zeros(27,3);
stdev = zeros(27,3);
wvrms = zeros(27,1);
depth = zeros(27,1);
e = zeros(3,3);avg = zeros(3,3);wrms = zeros(3,1);
tide = zeros(3,1);stdv = zeros(3,3);

%wvfiles are split by experiment, not by day
whichwv = [3 3 4 4 5 6 7 8 9];
ebb = [1 0 0 1 0 0 1 0 0];                                                  %experiments where the ebb tide was recorded
c = 1;cc = 1;
for i = 3:length(tkefiles)
    load([tkedir tkefiles{i}])
    load([veldir velfiles{i}])
    k = whichwv(c);load([wvdir wvfiles{k}])
    name{c} = regexprep(velfiles{i},'Vave.mat','');
    
    %some experiments lasted the whole tidal cycle. Crop to low-high or
    %high-low tide
    if i == 9
        start = datenum(2015,03,09,04,30,00);
    else
        start = Stat.vpro1.time(1);
    end
    if i == 8
        stop = datenum(2015,03,08,16,48,00);
    elseif i == 11
        stop = datenum(2015,03,14,09,20,00);
    else
        stop = Stat.vpro1.time(end);
    end
    %Use times to determine the number of samples to average over by
    %dividing the timeseries in thirds
    [~,~,~,hr,mi,~] = datevec(stop-start);
    time = (hr*60) + mi;                                                    %calculate minutes elapsed between start/stop time of experiment
    time = datenum(0,0,0,0,floor(time/3),0);
    if ebb(c) == 0
        timee = start:time:stop;timee(end) = stop;
    elseif ebb(c) == 1
        timee = fliplr(start:time:stop); timee(1) = stop;
    end
    for j = 1:length(timee)-1
        if ebb(c) == 0
            id = find(Stat.vpro1.time >= timee(j) & Stat.vpro1.time <= timee(j+1));
            id2 = find(wvstats.time >= timee(j) & wvstats.time <= timee(j+1));
        elseif ebb(c) == 1
            id = find(Stat.vpro1.time <= timee(j) & Stat.vpro1.time >= timee(j+1));
            id2 = find(wvstats.time <= timee(j) & wvstats.time >= timee(j+1));
        end
        tk(1) = nanmean((Stat.vpro1.beam1.E(5,id)+Stat.vpro1.beam3.E(5,id))./2);
        av(1) = nanmean(Avgs.vpro1.x(5,id));
        sd(1) = std((Stat.vpro1.beam1.E(5,id)+Stat.vpro1.beam3.E(5,id))./2);
        tk(2) = nanmean((Stat.vpro2.beam1.E(5,id)+Stat.vpro2.beam3.E(5,id))./2);
        av(2) = nanmean(Avgs.vpro2.x(5,id));
        sd(2) = std((Stat.vpro2.beam1.E(5,id)+Stat.vpro2.beam3.E(5,id))./2);
        tk(3) = nanmean((Stat.vpro3.beam1.E(5,id)+Stat.vpro3.beam3.E(5,id))./2);
        av(3) = nanmean(Avgs.vpro3.x(5,id));
        sd(3) = std((Stat.vpro3.beam1.E(5,id)+Stat.vpro3.beam3.E(5,id))./2);
        hrms = nanmean(wvstats.hrmsp(id2));
        dep = nanmean(wvstats.depth(id2));
        %save to variable
        e(j,:) = tk;
        avg(j,:) = av;
        wrms(j,:) = hrms;
        tide(j,:) = dep;
        stdv(j,:) = sd;
    end
    eps(cc:cc+2,:) = e;
    avgs(cc:cc+2,:) = abs(avg);
    wvrms(cc:cc+2,:) = wrms;
    depth(cc:cc+2,:) = tide;
    stdev(cc:cc+2,:) = stdv;
    c = c+1;
    cc = cc+3;
end
E = NaN(27,3);
E(1:3,:) = eps(1:3,:);E(4:6,:) = eps(4:6,:);E(7:9,:) = eps(13:15,:);
E(10:12,:) = eps(16:18,:);E(13:15,:) = eps(19:21,:);E(16:18,:) = eps(22:24,:);
E(19:21,:) = eps(7:9,:);E(22:24,:) = eps(10:12,:);E(25:27,:) = eps(25:27,:);
Avg = NaN(27,3);
Avg(1:3,:) = avgs(1:3,:);Avg(4:6,:) = avgs(4:6,:);Avg(7:9,:) = avgs(13:15,:);
Avg(10:12,:) = avgs(16:18,:);Avg(13:15,:) = avgs(19:21,:);Avg(16:18,:) = avgs(22:24,:);
Avg(19:21,:) = avgs(7:9,:);Avg(22:24,:) = avgs(10:12,:);Avg(25:27,:) = avgs(25:27,:);
Depth = NaN(27,1);
Depth(1:3,:) = depth(1:3,:);Depth(4:6,:) = depth(4:6,:);Depth(7:9,:) = depth(13:15,:);
Depth(10:12,:) = depth(16:18,:);Depth(13:15,:) = depth(19:21,:);Depth(16:18,:) = depth(22:24,:);
Depth(19:21,:) = depth(7:9,:);Depth(22:24,:) = depth(10:12,:);Depth(25:27,:) = depth(25:27,:);
Depth = repmat(Depth,1,3);
Wrms = NaN(27,1);
Wrms(1:3,:) = wvrms(1:3,:);Wrms(4:6,:) = wvrms(4:6,:);Wrms(7:9,:) = wvrms(13:15,:);
Wrms(10:12,:) = wvrms(16:18,:);Wrms(13:15,:) = wvrms(19:21,:);Wrms(16:18,:) = wvrms(22:24,:);
Wrms(19:21,:) = wvrms(7:9,:);Wrms(22:24,:) = wvrms(10:12,:);Wrms(25:27,:) = wvrms(25:27,:);
Wrms = repmat(Wrms,1,3);
Stdev = NaN(27,3);
Stdev(1:3,:) = stdev(1:3,:);Stdev(4:6,:) = stdev(4:6,:);Stdev(7:9,:) = stdev(13:15,:);
Stdev(10:12,:) = stdev(16:18,:);Stdev(13:15,:) = stdev(19:21,:);Stdev(16:18,:) = stdev(22:24,:);
Stdev(19:21,:) = stdev(7:9,:);Stdev(22:24,:) = stdev(10:12,:);Stdev(25:27,:) = stdev(25:27,:);
quad = cell(27,3);
quad(1:3,:) = name(1);quad(4:6,:) = name(2);quad(7:9,:) = name(5);
quad(10:12,:) = name(6);quad(13:15,:) = name(7);quad(16:18,:) = name(8);
quad(19:21,:) = name(3);quad(22:24,:) = name(4);quad(25:27,:) = name(9);
veldat.five.quad = quad;veldat.five.E = E;veldat.five.Estd = Stdev;
veldat.five.avg = Avg;veldat.five.depth = Depth;veldat.five.wrms = Wrms;

clearvars -except veldat vegdat