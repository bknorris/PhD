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

%Update: 23/06/2016: script updated to account for changes to the TKE, Avg
%and Wvstats files recreated mid June '16. Note that TKE is averaged every
%10s, while Avg is averaged every 30s.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%% Load the 2014 Data %%%%
veldir = 'd:\Projects\Mekong_F2014\DataAnalysis\Paper1\AveragedVelocities\';
wvdir = 'd:\Projects\Mekong_F2014\DataAnalysis\Paper1\Spectra\';
tkedir = 'D:\Projects\Mekong_F2014\DataAnalysis\Paper1\TKE\Vertical\';
velfiles = dir(veldir);velfiles = {velfiles.name};
tkefiles = dir(tkedir);tkefiles = {tkefiles.name};
wvfiles = dir(wvdir);wvfiles = {wvfiles.name};
%%%%extract TKE estimates from TKE files%%%%

%%%Extraction Settings%%%
horiz = [1 1 0 0 1 0 1 0];
ebb = [0 1 0 1 0 0 0 0];                                                    %experiments where the ebb tide was recorded
eps = zeros(24,4);
avgs = zeros(24,4);
ucubed = zeros(24,4);
ucur = zeros(24,4);
uwav = zeros(24,4);
umag = zeros(24,4);
uwstd = zeros(24,4);
stdev = zeros(24,4);
wvrms = zeros(24,4);
depth = zeros(24,4);
name = cell(8,1);
c = 1;cc = 1;                                                               %index counters
for i = 3:length(tkefiles)
    e = zeros(4,3);
    avg = zeros(4,3);
    ucub = zeros(4,3);
    Ucur = zeros(4,3);
    Uwav = zeros(4,3);
    Umagn = zeros(4,3);
    Uwstdv = zeros(4,3);
    wrms = zeros(4,3);
    tide = zeros(4,3);
    stdv = zeros(4,3);
    load([tkedir tkefiles{i}])
    load([veldir velfiles{i}])
    load([wvdir wvfiles{i}])
    name{c} = regexprep(velfiles{i},'Vave.mat','');

    %determine how many instruments are in the file
    fn = fieldnames(Avgs);
    [tp,~] = size(fn);
    if tp == 1
        g = 1;
    else
        g = 1:3;
    end
    for k = g                                                               %loop through instruments
        %%%Average TKE to the same length as Avgs
%         n = length(Avgs.(fn{k}).time);
%         beam1 = zeros(30,n);beam2 = zeros(30,n);beam3 = zeros(30,n);beam4 = zeros(30,n);
%         for kk = 1:n-1
%             td = find(Stat.(fn{k}).time >= Avgs.(fn{k}).time(kk) & Stat.(fn{k}).time <= Avgs.(fn{k}).time(kk+1));
%             beam1(1:30,kk) = nanmean(Stat.(fn{k}).beam1.E(:,td),2);
%             beam2(1:30,kk) = nanmean(Stat.(fn{k}).beam2.E(:,td),2);
%             beam3(1:30,kk) = nanmean(Stat.(fn{k}).beam3.E(:,td),2);
%             beam4(1:30,kk) = nanmean(Stat.(fn{k}).beam4.E(:,td),2);
%         end
        beam1 = Stat.(fn{k}).z1.E;beam2 = Stat.(fn{k}).z2.E;
        %%%Use times to determine the number of samples to average over by
        %dividing the timeseries into four parts (LL, ML, MH, HH)
        start = Stat.(fn{k}).time(1);stop = Stat.(fn{k}).time(end);
        [~,~,~,hr,mi,~] = datevec(stop-start);
        time = (hr*60) + mi;                                                %calculate minutes elapsed between start/stop time of experiment
        time = datenum(0,0,0,0,floor(time/4),0);
        if ebb(c) == 0
            timee = start:time:stop;timee(end) = stop;                      %if the data is from an ebb tide, flip the time record around so the first measurement is lower tidal elevation
        elseif ebb(c) == 1
            timee = fliplr(start:time:stop);timee(1) = stop;
        end
        %%%Average into subsections of the tide:
        for j = 1:length(timee)-1
            if ebb(c) == 0
                id = find(Avgs.(fn{k}).time >= timee(j) & Avgs.(fn{k}).time <= timee(j+1));
                id2 = find(wvstats.time >= timee(j) & wvstats.time <= timee(j+1));
            elseif ebb(c) == 1
                id = find(Avgs.(fn{k}).time <= timee(j) & Avgs.(fn{k}).time >= timee(j+1));
                id2 = find(wvstats.time <= timee(j) & wvstats.time >= timee(j+1));
            end
            if horiz(c) == 1                                                %if the instrument is horizontal, use bin 15 for TKE estimates
                tk = nanmean((beam1(15,id)+beam2(15,id))./2);
                uc = nanmean(Avgs.(fn{k}).ucubed(15,id));
                av = nanmean(Avgs.(fn{k}).x(15,id));
                Uc = nanmean(Avgs.(fn{k}).Uc(15,id));
                Uw = nanmean(Avgs.(fn{k}).Uw(15,id));
                Uwstd = nanstd(Avgs.(fn{k}).Uw(15,id));
                Umag = nanmean(Avgs.(fn{k}).Umag(15,id));
                sd = nanstd((beam1(15,id)+beam2(15,id))./2);
                hrms = nanmean(wvstats.hrmsp(id2));
                dep = nanmean(wvstats.depth(id2));
            elseif horiz(c) == 0
                tk = nanmean((beam1(5,id)+beam2(5,id))./2);
                uc = nanmean(Avgs.(fn{k}).ucubed(5,id));
                av = nanmean(Avgs.(fn{k}).x(5,id));
                Uc = nanmean(Avgs.(fn{k}).Uc(5,id));
                Uw = nanmean(Avgs.(fn{k}).Uw(5,id));
                Uwstd = nanstd(Avgs.(fn{k}).Uw(15,id));
                Umag = nanmean(Avgs.(fn{k}).Umag(5,id));
                sd = nanstd((beam1(5,id)+beam2(5,id))./2);
                hrms = nanmean(wvstats.hrmsp(id2));
                dep = nanmean(wvstats.depth(id2));
            end
            if j > 4                                                        %the FSS_1 experiment is so short that it has >4 timesteps
                continue
            else 
                e(j,k) = tk;
                avg(j,k) = av;
                ucub(j,k) = uc;
                wrms(j,k) = hrms;
                Ucur(j,k) = Uc;
                Uwav(j,k) = Uw;
                Uwstdv(j,k) = Uwstd;
                Umagn(j,k) = Umag;
                tide(j,k) = dep;
                stdv(j,k) = sd;
            end
        end
    end
    
    eps(cc:cc+2,:) = e';
    avgs(cc:cc+2,:) = avg';
    ucubed(cc:cc+2,:) = ucub';
    ucur(cc:cc+2,:) = Ucur';
    uwav(cc:cc+2,:) = Uwav';
    uwstd(cc:cc+2,:) = Uwstdv';
    umag(cc:cc+2,:) = Umagn';
    wvrms(cc:cc+2,:) = wrms';
    depth(cc:cc+2,:) = tide';
    stdev(cc:cc+2,:) = stdv';
    c = c+1;
    cc = cc+3;
end
%rearrange the data in the order of deployment (Q2-Q7).
E = NaN(24,4);
E(1:3,:) = eps(22:24,:);E(4:6,:) = eps(16:18,:);E(7:9,:) = eps(19:21,:);
E(10:12,:) = eps(1:3,:);E(13:15,:) = eps(13:15,:);E(16:18,:) = eps(7:9,:);
E(19:21,:) = eps(10:12,:);E(22:24,:) = eps(4:6,:);
Avg = NaN(24,4);
Avg(1:3,:) = avgs(22:24,:);Avg(4:6,:) = avgs(16:18,:);Avg(7:9,:) = avgs(19:21,:);
Avg(10:12,:) = avgs(1:3,:);Avg(13:15,:) = avgs(13:15,:);Avg(16:18,:) = avgs(7:9,:);
Avg(19:21,:) = avgs(10:12,:);Avg(22:24,:) = avgs(4:6,:);
UC = NaN(24,4);
UC(1:3,:) = ucubed(22:24,:);UC(4:6,:) = ucubed(16:18,:);UC(7:9,:) = ucubed(19:21,:);
UC(10:12,:) = ucubed(1:3,:);UC(13:15,:) = ucubed(13:15,:);UC(16:18,:) = ucubed(7:9,:);
UC(19:21,:) = ucubed(10:12,:);UC(22:24,:) = ucubed(4:6,:);
UCUR = NaN(24,4);
UCUR(1:3,:) = ucur(22:24,:);UCUR(4:6,:) = ucur(16:18,:);UCUR(7:9,:) = ucur(19:21,:);
UCUR(10:12,:) = ucur(1:3,:);UCUR(13:15,:) = ucur(13:15,:);UCUR(16:18,:) = ucur(7:9,:);
UCUR(19:21,:) = ucur(10:12,:);UCUR(22:24,:) = ucur(4:6,:);
UWAV = NaN(24,4);
UWAV(1:3,:) = uwav(22:24,:);UWAV(4:6,:) = uwav(16:18,:);UWAV(7:9,:) = uwav(19:21,:);
UWAV(10:12,:) = uwav(1:3,:);UWAV(13:15,:) = uwav(13:15,:);UWAV(16:18,:) = uwav(7:9,:);
UWAV(19:21,:) = uwav(10:12,:);UWAV(22:24,:) = uwav(4:6,:);
UMAG = NaN(24,4);
UMAG(1:3,:) = umag(22:24,:);UMAG(4:6,:) = umag(16:18,:);UMAG(7:9,:) = umag(19:21,:);
UMAG(10:12,:) = umag(1:3,:);UMAG(13:15,:) = umag(13:15,:);UMAG(16:18,:) = umag(7:9,:);
UMAG(19:21,:) = umag(10:12,:);UMAG(22:24,:) = umag(4:6,:);
UWSTD = NaN(24,4);
UWSTD(1:3,:) = uwstd(22:24,:);UWSTD(4:6,:) = uwstd(16:18,:);UWSTD(7:9,:) = uwstd(19:21,:);
UWSTD(10:12,:) = uwstd(1:3,:);UWSTD(13:15,:) = uwstd(13:15,:);UWSTD(16:18,:) = uwstd(7:9,:);
UWSTD(19:21,:) = uwstd(10:12,:);UWSTD(22:24,:) = uwstd(4:6,:);
Depth = NaN(24,4);
Depth(1:3,:) = depth(22:24,:);Depth(4,:) = depth(17,:);Depth(5:6,:) = NaN;
Depth(7:9,:) = depth(19:21,:);Depth(10:12,:) = depth(1:3,:);
Depth(13:15,:) = depth(13:15,:);Depth(16:18,:) = depth(7:9,:);
Depth(19:21,:) = depth(10:12,:);Depth(22:24,:) = depth(4:6,:);
Wrms = NaN(24,4);
Wrms(1:3,:) = wvrms(22:24,:);Wrms(4,:) = wvrms(17,:);Wrms(5:6,:) = NaN;
Wrms(7:9,:) = wvrms(19:21,:);Wrms(10:12,:) = wvrms(1:3,:);
Wrms(13:15,:) = wvrms(13:15,:);Wrms(16:18,:) = wvrms(7:9,:);
Wrms(19:21,:) = wvrms(10:12,:);Wrms(22:24,:) = wvrms(4:6,:);
Stdev = zeros(24,4);
Stdev(1:3,:) = stdev(22:24,:);Stdev(4,:) = stdev(17,:);Stdev(5:6,:) = NaN;
Stdev(7:9,:) = stdev(19:21,:);Stdev(10:12,:) = stdev(1:3,:);
Stdev(13:15,:) = stdev(13:15,:);Stdev(16:18,:) = stdev(7:9,:);
Stdev(19:21,:) = stdev(10:12,:);Stdev(22:24,:) = stdev(4:6,:);
quad = cell(24,4);
quad(1:3,:) = name(8);quad(4:6,:) = name(6);
quad(7:9,:) = name(7);quad(10:12,:) = name(1);
quad(13:15,:) = name(5);quad(16:18,:) = name(3);
quad(19:21,:) = name(4);quad(22:24,:) = name(2);
veldat.four.quad = quad;veldat.four.E = E;veldat.four.Estd = Stdev;
veldat.four.avg = Avg;veldat.four.ucubed = UC;veldat.four.uc = UCUR;
veldat.four.uw = UWAV;veldat.four.umag = UMAG;veldat.four.uwstd = UWSTD;
veldat.four.depth = Depth;
veldat.four.wrms = Wrms;

clearvars -except veldat vegdat

%%%% Load the 2015 Data %%%%
veldir = 'd:\Projects\Mekong_W2015\DataAnalysis\Paper1\AveragedVelocities\Original\';
tkedir = 'd:\Projects\Mekong_W2015\DataAnalysis\Paper1\TKE\Vertical\';
wvdir = 'd:\Projects\Mekong_W2015\DataAnalysis\Paper1\Spectra\';
velfiles = dir(veldir);velfiles = {velfiles.name};
tkefiles = dir(tkedir);tkefiles = {tkefiles.name};
wvfiles = dir(wvdir);wvfiles = {wvfiles.name};

eps = zeros(30,4);
avgs = zeros(30,4);
ucubed = zeros(30,4);
ucur = zeros(30,4);
uwav = zeros(30,4);
uwstd = zeros(30,4);
umag = zeros(30,4);
stdev = zeros(30,4);
wvrms = zeros(30,4);
depth = zeros(30,4);
e = zeros(4,3);avg = zeros(4,3);wrms = zeros(4,3);
tide = zeros(4,3);stdv = zeros(4,3);

%wvfiles are split by experiment, not by day
whichwv = [3 3 6 6 8 9 10 10 10 11 12 12 12 12 12 12];
ebb = [1 0 0 1 0 0 1 1 1 0 1 1 1 0 0 0];                                                  %experiments where the ebb tide was recorded
c = 1;cc = 1;
for i = 3:length(tkefiles)
    e = zeros(4,3);
    avg = zeros(4,3);
    ucub = zeros(4,3);
    Ucur = zeros(4,3);
    Uwav = zeros(4,3);
    Uwstdv = zeros(4,3);
    Umagn = zeros(4,3);
    wrms = zeros(4,3);
    tide = zeros(4,3);
    stdv = zeros(4,3);
    load([tkedir tkefiles{i}])
    load([veldir velfiles{i}])
    wf = whichwv(c);load([wvdir wvfiles{wf}])
    name{c} = regexprep(velfiles{i},'Vave.mat','');
    
    %determine how many instruments are in the file
    fn = fieldnames(Avgs);
    [tp,~] = size(fn);
    if tp == 1
        g = 1;
    else
        g = 1:3;
    end
    for k = g                                                               %loop through instruments
        %%%Average TKE to the same length as Avgs
%         n = length(Avgs.(fn{k}).time);
%         beam1 = zeros(30,n);beam2 = zeros(30,n);beam3 = zeros(30,n);beam4 = zeros(30,n);
%         for kk = 1:n-1
%             td = find(Stat.(fn{k}).time >= Avgs.(fn{k}).time(kk) & Stat.(fn{k}).time <= Avgs.(fn{k}).time(kk+1));
%             beam1(1:30,kk) = nanmean(Stat.(fn{k}).beam1.E(:,td),2);
%             beam2(1:30,kk) = nanmean(Stat.(fn{k}).beam2.E(:,td),2);
%             beam3(1:30,kk) = nanmean(Stat.(fn{k}).beam3.E(:,td),2);
%             beam4(1:30,kk) = nanmean(Stat.(fn{k}).beam4.E(:,td),2);
%         end
        beam1 = Stat.(fn{k}).z1.E;beam2 = Stat.(fn{k}).z2.E;
        %some experiments lasted the whole tidal cycle. Crop to low-high or
        %high-low tide
        if i >= 9 && i <= 11
            start = datenum(2015,03,09,04,30,00);
        else
            start = Stat.(fn{k}).time(1);
        end
        if i == 5
            stop = datenum(2015,03,11,16,00,00);
        elseif i == 8
            stop = datenum(2015,03,08,16,48,00);
        elseif i == 15
            stop = datenum(2015,03,13,22,00,00);
        elseif i >= 16 && i <= 18
            stop = datenum(2015,03,14,09,20,00);
        else
            stop = Stat.(fn{k}).time(end);
        end
        %%%Use times to determine the number of samples to average over by
        %dividing the timeseries into four parts (LL, ML, MH, HH)
        [~,~,~,hr,mi,~] = datevec(stop-start);
        time = (hr*60) + mi;                                                %calculate minutes elapsed between start/stop time of experiment
        time = datenum(0,0,0,0,floor(time/4),0);
        if ebb(c) == 0
            timee = start:time:stop;timee(end) = stop;                      %if the data is from an ebb tide, flip the time record around so the first measurement is lower tidal elevation
        elseif ebb(c) == 1
            timee = fliplr(start:time:stop);timee(1) = stop;
        end
        %%%Average into subsections of the tide:
        for j = 1:length(timee)-1
            if ebb(c) == 0
                id = find(Avgs.(fn{k}).time >= timee(j) & Avgs.(fn{k}).time <= timee(j+1));
                id2 = find(wvstats.time >= timee(j) & wvstats.time <= timee(j+1));
            elseif ebb(c) == 1
                id = find(Avgs.(fn{k}).time <= timee(j) & Avgs.(fn{k}).time >= timee(j+1));
                id2 = find(wvstats.time <= timee(j) & wvstats.time >= timee(j+1));
            end
            tk = nanmean((beam1(5,id)+beam2(5,id))./2);
            uc = nanmean(Avgs.(fn{k}).ucubed(5,id));
            av = nanmean(Avgs.(fn{k}).x(5,id));
            Uc = nanmean(Avgs.(fn{k}).Uc(15,id));
            Uw = nanmean(Avgs.(fn{k}).Uw(15,id));
            Uwstd = nanstd(Avgs.(fn{k}).Uw(15,id));
            Umag = nanmean(Avgs.(fn{k}).Umag(15,id));
            sd = nanstd((beam1(5,id)+beam2(5,id))./2);
            hrms = nanmean(wvstats.hrmsp(id2));
            dep = nanmean(wvstats.depth(id2));
            if j > 4                                                        %the FSS_1 experiment is so short that it has >4 timesteps
                continue
            else
                e(j,k) = tk;
                avg(j,k) = av;
                ucub(j,k) = uc;
                wrms(j,k) = hrms;
                Ucur(j,k) = Uc;
                Uwav(j,k) = Uw;
                Uwstdv(j,k) = Uwstd;
                Umagn(j,k) = Umag;
                tide(j,k) = dep;
                stdv(j,k) = sd;
            end
        end
    end
    if length(g) > 1
        eps(cc:cc+2,:) = e';
        avgs(cc:cc+2,:) = avg';
        ucubed(cc:cc+2,:) = ucub';
        ucur(cc:cc+2,:) = Ucur';
        uwav(cc:cc+2,:) = Uwav';
        uwstd(cc:cc+2,:) = Uwstdv';
        umag(cc:cc+2,:) = Umagn';
        wvrms(cc:cc+2,:) = wrms';
        depth(cc:cc+2,:) = tide';
        stdev(cc:cc+2,:) = stdv';
        c = c+1;
        cc = cc+3;
    else
        eps(cc,:) = e(:,1)';
        avgs(cc,:) = avg(:,1)';
        ucubed(cc,:) = ucub(:,1)';
        ucur(cc,:) = Ucur(:,1)';
        uwav(cc,:) = Uwav(:,1)';
        uwstd(cc,:) = Uwstdv(:,1)';
        umag(cc,:) = Umagn(:,1)';
        wvrms(cc,:) = wrms(:,1)';
        depth(cc,:) = tide(:,1)';
        stdev(cc,:) = stdv(:,1)';
        c = c+1;
        cc = cc+1;
    end
end
E = NaN(30,4);
E(1:3,:) = eps(1:3,:);E(4:6,:) = eps(4:6,:);E(7:9,:) = eps(13:15,:);
E(10:12,:) = eps(16:18,:);E(13,:) = eps(19,:);E(14,:) = eps(20,:);
E(15,:) = eps(21,:);E(16:18,:) = eps(22:24,:);E(19:21,:) = eps(7:9,:);
E(22:24,:) = eps(10:12,:);E(25,:) = eps(25,:);E(26,:) = eps(26,:);
E(27,:) = eps(27,:);E(28,:) = eps(28,:);E(29,:) = eps(29,:);E(30,:) = eps(30,:);
Avg = NaN(30,4);
Avg(1:3,:) = avgs(1:3,:);Avg(4:6,:) = avgs(4:6,:);Avg(7:9,:) = avgs(13:15,:);
Avg(10:12,:) = avgs(16:18,:);Avg(13,:) = avgs(19,:);Avg(14,:) = avgs(20,:);
Avg(15,:) = avgs(21,:);Avg(16:18,:) = avgs(22:24,:);Avg(19:21,:) = avgs(7:9,:);
Avg(22:24,:) = avgs(10:12,:);Avg(25,:) = avgs(25,:);Avg(26,:) = avgs(26,:);
Avg(27,:) = avgs(27,:);Avg(28,:) = avgs(28,:);Avg(29,:) = avgs(29,:);Avg(30,:) = avgs(30,:);
UC = NaN(30,4);
UC(1:3,:) = ucubed(1:3,:);UC(4:6,:) = ucubed(4:6,:);UC(7:9,:) = ucubed(13:15,:);
UC(10:12,:) = ucubed(16:18,:);UC(13,:) = ucubed(19,:);UC(14,:) = ucubed(20,:);
UC(15,:) = ucubed(21,:);UC(16:18,:) = ucubed(22:24,:);UC(19:21,:) = ucubed(7:9,:);
UC(22:24,:) = ucubed(10:12,:);UC(25,:) = ucubed(25,:);UC(26,:) = ucubed(26,:);
UC(27,:) = ucubed(27,:);UC(28,:) = ucubed(28,:);UC(29,:) = ucubed(29,:);UC(30,:) = ucubed(30,:);
UCUR = NaN(30,4);
UCUR(1:3,:) = ucur(1:3,:);UCUR(4:6,:) = ucur(4:6,:);UCUR(7:9,:) = ucur(13:15,:);
UCUR(10:12,:) = ucur(16:18,:);UCUR(13,:) = ucur(19,:);UCUR(14,:) = ucur(20,:);
UCUR(15,:) = ucur(21,:);UCUR(16:18,:) = ucur(22:24,:);UCUR(19:21,:) = ucur(7:9,:);
UCUR(22:24,:) = ucur(10:12,:);UCUR(25,:) = ucur(25,:);UCUR(26,:) = ucur(26,:);
UCUR(27,:) = ucur(27,:);UCUR(28,:) = ucur(28,:);UCUR(29,:) = ucur(29,:);UCUR(30,:) = ucur(30,:);
UWAV = NaN(30,4);
UWAV(1:3,:) = uwav(1:3,:);UWAV(4:6,:) = uwav(4:6,:);UWAV(7:9,:) = uwav(13:15,:);
UWAV(10:12,:) = uwav(16:18,:);UWAV(13,:) = uwav(19,:);UWAV(14,:) = uwav(20,:);
UWAV(15,:) = uwav(21,:);UWAV(16:18,:) = uwav(22:24,:);UWAV(19:21,:) = uwav(7:9,:);
UWAV(22:24,:) = uwav(10:12,:);UWAV(25,:) = uwav(25,:);UWAV(26,:) = uwav(26,:);
UWAV(27,:) = uwav(27,:);UWAV(28,:) = uwav(28,:);UWAV(29,:) = uwav(29,:);UWAV(30,:) = uwav(30,:);
UMAG = NaN(30,4);
UMAG(1:3,:) = umag(1:3,:);UMAG(4:6,:) = umag(4:6,:);UMAG(7:9,:) = umag(13:15,:);
UMAG(10:12,:) = umag(16:18,:);UMAG(13,:) = umag(19,:);UMAG(14,:) = umag(20,:);
UMAG(15,:) = umag(21,:);UMAG(16:18,:) = umag(22:24,:);UMAG(19:21,:) = umag(7:9,:);
UMAG(22:24,:) = umag(10:12,:);UMAG(25,:) = umag(25,:);UMAG(26,:) = umag(26,:);
UMAG(27,:) = umag(27,:);UMAG(28,:) = umag(28,:);UMAG(29,:) = umag(29,:);UMAG(30,:) = umag(30,:);
UWSTD = NaN(30,4);
UWSTD(1:3,:) = uwstd(1:3,:);UWSTD(4:6,:) = uwstd(4:6,:);UWSTD(7:9,:) = uwstd(13:15,:);
UWSTD(10:12,:) = uwstd(16:18,:);UWSTD(13,:) = uwstd(19,:);UWSTD(14,:) = uwstd(20,:);
UWSTD(15,:) = uwstd(21,:);UWSTD(16:18,:) = uwstd(22:24,:);UWSTD(19:21,:) = uwstd(7:9,:);
UWSTD(22:24,:) = uwstd(10:12,:);UWSTD(25,:) = uwstd(25,:);UWSTD(26,:) = uwstd(26,:);
UWSTD(27,:) = uwstd(27,:);UWSTD(28,:) = uwstd(28,:);UWSTD(29,:) = uwstd(29,:);UWSTD(30,:) = uwstd(30,:);
Depth = NaN(30,4);
Depth(1:3,:) = depth(1:3,:);Depth(4:6,:) = depth(4:6,:);Depth(7:9,:) = depth(13:15,:);
Depth(10:12,:) = depth(16:18,:);Depth(13,:) = depth(19,:);Depth(14,:) = depth(20,:);
Depth(15,:) = depth(21,:);Depth(16:18,:) = depth(22:24,:);Depth(19:21,:) = depth(7:9,:);
Depth(22:24,:) = depth(10:12,:);Depth(25,:) = depth(25,:);Depth(26,:) = depth(26,:);
Depth(27,:) = depth(27,:);Depth(28,:) = depth(28,:);Depth(29,:) = depth(29,:);Depth(30,:) = depth(30,:);
Wrms = NaN(30,4);
Wrms(1:3,:) = wvrms(1:3,:);Wrms(4:6,:) = wvrms(4:6,:);Wrms(7:9,:) = wvrms(13:15,:);
Wrms(10:12,:) = wvrms(16:18,:);Wrms(13,:) = wvrms(19,:);Wrms(14,:) = wvrms(20,:);
Wrms(15,:) = wvrms(21,:);Wrms(16:18,:) = wvrms(22:24,:);Wrms(19:21,:) = wvrms(7:9,:);
Wrms(22:24,:) = wvrms(10:12,:);Wrms(25,:) = wvrms(25,:);Wrms(26,:) = wvrms(26,:);
Wrms(27,:) = wvrms(27,:);Wrms(28,:) = wvrms(28,:);Wrms(29,:) = wvrms(29,:);Wrms(30,:) = wvrms(30,:);
Stdev = NaN(30,4);
Stdev(1:3,:) = stdev(1:3,:);Stdev(4:6,:) = stdev(4:6,:);Stdev(7:9,:) = stdev(13:15,:);
Stdev(10:12,:) = stdev(16:18,:);Stdev(13,:) = stdev(19,:);Stdev(14,:) = stdev(20,:);
Stdev(15,:) = stdev(21,:);Stdev(16:18,:) = stdev(22:24,:);Stdev(19:21,:) = stdev(7:9,:);
Stdev(22:24,:) = stdev(10:12,:);Stdev(25,:) = stdev(25,:);Stdev(26,:) = stdev(26,:);
Stdev(27,:) = stdev(27,:);Stdev(28,:) = stdev(28,:);Stdev(29,:) = stdev(29,:);Stdev(30,:) = stdev(30,:);
quad = cell(30,4);
quad(1:3,:) = name(1);quad(4:6,:) = name(2);quad(7:9,:) = name(5);
quad(10:12,:) = name(6);quad(13,:) = name(7);quad(14,:) = name(8);
quad(15,:) = name(9);quad(16:18,:) = name(10);quad(19:21,:) = name(3);
quad(22:24,:) = name(4);quad(25,:) = name(11);quad(26,:) = name(12);
quad(27,:) = name(13);quad(28,:) = name(14);quad(29,:) = name(15);quad(30,:) = name(16);
veldat.five.quad = quad;veldat.five.E = E;veldat.five.Estd = Stdev;
veldat.five.avg = Avg;veldat.five.ucubed = UC;veldat.five.uc = UCUR;
veldat.five.uw = UWAV;veldat.five.umag = UMAG;veldat.five.uwstd = uwstd;
veldat.five.depth = Depth;
veldat.five.wrms = Wrms;

clearvars -except veldat vegdat
save(['d:\Projects\Mekong_W2015\DataAnalysis\Paper1\Environment\' 'VelWvTKEdata'],'veldat','-v7.3')
