%Average upper-canopy ADVs from the HTA. Note that because the pressure
%records for these ADVs are screwy (HTA1-2), I'll also load the other ADV,
%sync the time series up, and use those for depth calculations.
clear
datdir = 'd:\Projects\Mekong_W2015\Data\Vector\FSS\';
vc1 = {'VC01_070315.mat';'VC01_080315.mat';''};
vc5 = {'V5109_070315.mat';'V5109_080315.mat';'V5109_100315.mat'};
%Crop times
cstart = {'07-Mar-2015 13:36:32';'08-Mar-2015 14:14:22';'10-Mar-2015 14:45:06'};
cstop = {'07-Mar-2015 17:08:02';'08-Mar-2015 19:04:22';'10-Mar-2015 16:39:36'};
Vave = struct();
dates = {'day1';'day2';'day3'};
heading = [340 340 340]; %headings are in order! (HTA1 - 3)
for i = 1:3
    disp(['Running ' dates{i}])
    %Load the data
    if ~isempty(vc1{i})
        load([datdir vc1{i}]) %VC01
        t1 = ADV.datetime;
        vid = find(t1 >= datenum(cstart{i}) & t1 <= datenum(cstop{i}));        
        U = ADV.U(vid);
        V = ADV.V(vid);
        W = ADV.W(vid);
        t1 = t1(vid);
        load([datdir vc5{i}]) %V5109
        t2 = ADV.datetime;
        vid = find(t2 >= datenum(cstart{i}) & t2 <= datenum(cstop{i}));        
        P = ADV.Pres(vid)+ADV.Metadata.pressure_sensor_height/1000;
        t2 = t2(vid);
        ll = (sin(ADV.Metadata.inst_lat/57.29578))^2;
        fs = ADV.Metadata.instmeta.samprate;
        clear ADV
        %%%
    else
        load([datdir vc5{i}]) %V5109
        t1 = ADV.datetime;
        vid = find(t1 >= datenum(cstart{i}) & t1 <= datenum(cstop{i}));        
        U = ADV.U(vid);
        V = ADV.V(vid);
        W = ADV.W(vid);
        t1 = t1(vid);
        P = ADV.Pres(vid)+ADV.Metadata.pressure_sensor_height/1000;
        x = (sin(ADV.Metadata.inst_lat/57.29578))^2;
        fs = ADV.Metadata.instmeta.samprate;
        clear ADV
        %%%
    end
    %Averaging
    step = 120; %seconds
    win = 600; %seconds
    avt = fs*step;
    nwin = fs*win;
    %%%Data Analysis
    nsamp = length(vid);
    ind = [1 avt:avt:nsamp];
    for ii = 1:length(ind)
        if abs(nsamp-ind(ii)) < nwin  %skip the last few indexes approaching the end of the t-s
            continue
        else
            idx = ind(ii):ind(ii)+nwin-1;
        end
        u = U(idx);
        v = V(idx);
        w = W(idx);
        p = P(idx);
        time = t1(ind(ii));
        %Calculate depth
        g = zeros(length(p),1);h = zeros(length(p),1);
        for j = 1:length(p)
            g(j,:) = 9.780318*(1.0+(5.2788E-3+2.36E-5*ll)*ll)+1.092E-6*p(j,:);
            h(j,:) = ((((-1.82E-15*p(j,:)+2.279E-10)*p(j,:)-2.2512E-5)*p(j,:)+9.72659)*p(j,:))/g(j,:);
        end
        h = mean(h);
        g = mean(g);
        %rotate vels to along current and across current
        th = heading(i)*pi/180;
        R = [cos(th) -sin(th); sin(th) cos(th)];
        vxy = [v u];rxy = zeros(size(vxy)); %v is north, u is east
        for jj = 1:length(vxy),rxy(jj,:) = vxy(jj,:)*R;end
        x = vxy(:,1);y = vxy(:,2);x(isnan(x)) = [];y(isnan(y)) = [];
        %Calculate Uw and Uc, Luhar et al. 2013
        Ec = (1/length(y))*sum(y);Nc = (1/length(x))*sum(x);
        Ewrms = sqrt((1/length(y))*sum((y-Ec).^2));
        Nwrms = sqrt((1/length(x))*sum((x-Nc).^2));
        Uc = sqrt(Ec^2+Nc^2);Uwrms = sqrt(Ewrms^2+Nwrms^2);
        Uw = sqrt(2)*Uwrms;
        
            %save data to structure
        Vave.(dates{i}).time(ii) = time;
        Vave.(dates{i}).u(ii) = nanmean(x); %across-shore
        Vave.(dates{i}).v(ii) = nanmean(y); %along-shore
        Vave.(dates{i}).w(ii) = nanmean(w);
        Vave.(dates{i}).Uw(ii) = Uw;
        Vave.(dates{i}).Uc(ii) = Uc;
        Vave.(dates{i}).g(ii) = g;
        Vave.(dates{i}).h(ii) = h;
    end
end
savedatdir = 'd:\Projects\Mekong_W2015\DataAnalysis\Paper2\AveragedVelocities\';
save([savedatdir 'HTA_ADVave'],'Vave','-v7.3')