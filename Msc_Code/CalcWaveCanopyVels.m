%Plot wave and velocities for measurements at the bed and above the canopy,
%calculated from wvstats using the Vectrino and a collocated pressure
%sensor (either the ADV or the RBR Duet): compare measurements from HTA day
%1 for VPs 2 and 3 with ADV5109, then compare the synoptic measurements from
%the VTA, VPs 1 and 3.
clear
%%%%load day 1
load('C:\Users\bkn5\Projects\Mekong_W2015\DataAnalysis\Paper1\HTAday1Vels.mat')
%prep data for wvstats
heading = 20;
u = HTA.vpro2.x(:,1);
v = HTA.vpro2.y(:,1);
%rotate VP2 to the cross-shore direction
rot = heading*pi/180;
u2 = u.*(ones(size(u))*cos(rot)) + ...
    v.*(ones(size(v))*sin(rot));
v2 = -u.*(ones(size(u))*sin(rot)) + ...
    v.*(ones(size(v))*cos(rot));
%downsample to 2Hz
u2 = downsample(u2,25);
v2 = downsample(v2,25);
u = HTA.vpro3.x(:,1);
v = HTA.vpro3.y(:,1);
%rotate VP3 to the cross-shore direction
rot = heading*pi/180;
u3 = u.*(ones(size(u))*cos(rot)) + ...
    v.*(ones(size(v))*sin(rot));
v3 = -u.*(ones(size(u))*sin(rot)) + ...
    v.*(ones(size(v))*cos(rot));
%downsample to 2Hz
u3 = downsample(u3,25);bb = length(u3);
v3 = downsample(v3,25);
Uvp = (u2(1:bb)+u3)./2;Vvp = (v2(1:bb)+v3)./2;
clear HTA
%get pressure signal from collocated ADV
load('C:\Users\bkn5\Projects\Mekong_W2015\DataAnalysis\Paper1\HTAtimes.mat')
load('C:\Users\bkn5\Projects\Mekong_W2015\Data\Vector\FSS\V5109_070315.mat')
start = HTA.times.t1;stop = HTA.times.e1;
%crop data fields to the specified start/stop times of the experiment
ind = find(ADV.datetime >= start & ADV.datetime <=stop);
p = downsample(ADV.Pres(ind),16);
time = downsample(ADV.datetime(ind),16);
p = p(1:bb);time = time(1:bb);
x = (sin(ADV.Metadata.inst_lat/57.29578))^2;

%Iterate through time intervals
%wavestats settings
zp = 0.240; %height of ADV pressure sensor
zuv = 0.0211; %height of velocity measurement
fs = 2;
avt = fs*2*60; %# samples in 2 min
idx = avt:avt:length(p);
idxx = [1 idx];
it = 1;
while it < length(idxx)
    indx = idxx(it):idxx(it+1);
    P = p(indx);
    U = Uvp(indx);
    V = Vvp(indx);
    ts = time(idxx(it)); %save timestamp of the spectra interval
    if rem(length(P), 2) == 1
        n = length(P)-1; %force inputs to be even
        P = P(1:end-1);
        U = U(1:end-1);
        V = V(1:end-1);
    else
        n = length(P);
    end
    U(isnan(U)) = 0;
    V(isnan(V)) = 0;
    
    %Spectra settings
    nfft = floor(0.25*n);
    window = hanning(n,'periodic');
    min_f = 0.01; % mininum frequency, below which no correction is applied
    max_f = 1.5; % maximum frequency, above which no correction is applied
    
    g = 9.81;
    rho = 1.009629812357770e+03;
    for i = 1:length(p)
        g(i,:) = 9.780318*(1.0+(5.2788E-3+2.36E-5*x)*x)+1.092E-6*p(i,:);
        h(i,:) = ((((-1.82E-15*p(i,:)+2.279E-10)*p(i,:)-2.2512E-5)*p(i,:)+9.72659)*p(i,:))/g(i,:);
    end
    h = mean(h);
    ws = puv(P,U,V,h,zp,zuv,fs,nfft,rho,window,min_f,max_f);
    wvstats.hta.time(it,1) = time(idxx(it));
    wvstats.hta.wavev(it,1) = ws.ubr;
    wvstats.hta.azm(it,1) = ws.azr;
    it = it+1;
end

%%%%load VTA VP1
load('C:\Users\bkn5\Projects\Mekong_W2015\DataAnalysis\Paper1\VTAvelocities.mat')
%prep data for wvstats
heading = 100;
u = VTA.vpro1.x(:,7);
v = VTA.vpro1.y(:,7);
%rotate VP2 to the cross-shore direction
rot = heading*pi/180;
u = u.*(ones(size(u))*cos(rot)) + ...
    v.*(ones(size(v))*sin(rot));
v = -u.*(ones(size(u))*sin(rot)) + ...
    v.*(ones(size(v))*cos(rot));
%downsample to 2Hz
Uvp = downsample(u,25);bb = length(Uvp);
Vvp = downsample(v,25);
clear VTA

%get pressure signal from collocated Duet
load('C:\Users\bkn5\Projects\Mekong_W2015\DataAnalysis\Paper1\VTAtimes.mat')
load('C:\Users\bkn5\Projects\Mekong_W2015\Data\RBR\Duet\DPS2\Duet_140315.mat')
start = VTA.times.t2;stop = VTA.times.e2;
%crop data fields to the specified start/stop times of the experiment
ind = find(RBR.Datetime >= start & RBR.Datetime <=stop);
p = downsample(RBR.SeaPres(ind),8);
time = downsample(RBR.Datetime(ind),8);
p = p(1:bb);time = time(1:bb);
x = (sin(RBR.Metadata.latitude/57.29578))^2;

%Iterate through time intervals
%wavestats settings
zp = 0.52; %height of RBR pressure sensor
zuv = 0.0242; %height of velocity measurement
fs = 1;
avt = fs*2*60; %# samples in 2 min
idx = avt:avt:length(p);
idxx = [1 idx];
it = 1;
while it < length(idxx)
    indx = idxx(it):idxx(it+1);
    P = p(indx);
    U = Uvp(indx);
    V = Vvp(indx);
    ts = time(idxx(it)); %save timestamp of the spectra interval
    if rem(length(P), 2) == 1
        n = length(P)-1; %force inputs to be even
        P = P(1:end-1);
        U = U(1:end-1);
        V = V(1:end-1);
    else
        n = length(P);
    end
    U(isnan(U)) = 0;
    V(isnan(V)) = 0;
    
    %Spectra settings
    nfft = floor(0.25*n);
    window = hanning(n,'periodic');
    min_f = 0.01; % mininum frequency, below which no correction is applied
    max_f = 1.5; % maximum frequency, above which no correction is applied
    
    g = 9.81;
    rho = 1.009629812357770e+03;
    for i = 1:length(p)
        g(i,:) = 9.780318*(1.0+(5.2788E-3+2.36E-5*x)*x)+1.092E-6*p(i,:);
        h(i,:) = ((((-1.82E-15*p(i,:)+2.279E-10)*p(i,:)-2.2512E-5)*p(i,:)+9.72659)*p(i,:))/g(i,:);
    end
    h = mean(h);
    ws = puv(P,U,V,h,zp,zuv,fs,nfft,rho,window,min_f,max_f);
    wvstats.vta.vp1.time(it,1) = time(idxx(it));
    wvstats.vta.vp1.wavev(it,1) = ws.ubr;
    wvstats.vta.vp1.azm(it,1) = ws.azr;
    it = it+1;
end

%%%%load VTA VP3
load('C:\Users\bkn5\Projects\Mekong_W2015\DataAnalysis\Paper1\VTAvelocities.mat')
%prep data for wvstats
heading = 100;
u = VTA.vpro3.x(:,30);
v = VTA.vpro3.y(:,30);
%rotate VP2 to the cross-shore direction
rot = heading*pi/180;
u = u.*(ones(size(u))*cos(rot)) + ...
    v.*(ones(size(v))*sin(rot));
v = -u.*(ones(size(u))*sin(rot)) + ...
    v.*(ones(size(v))*cos(rot));
%downsample to 2Hz
Uvp = downsample(u,25);bb = length(Uvp);
Vvp = downsample(v,25);
clear VTA

%get pressure signal from collocated Duet
p = downsample(RBR.SeaPres(ind),8);
time = downsample(RBR.Datetime(ind),8);
p = p(1:bb);time = time(1:bb);
x = (sin(RBR.Metadata.latitude/57.29578))^2;

%Iterate through time intervals
%wavestats settings
zp = 0.52; %height of RBR pressure sensor
zuv = 0.73; %height of velocity measurement
fs = 1;
avt = fs*2*60; %# samples in 2 min
idx = avt:avt:length(p);
idxx = [1 idx];
it = 1;
while it < length(idxx)
    indx = idxx(it):idxx(it+1);
    P = p(indx);
    U = Uvp(indx);
    V = Vvp(indx);
    ts = time(idxx(it)); %save timestamp of the spectra interval
    if rem(length(P), 2) == 1
        n = length(P)-1; %force inputs to be even
        P = P(1:end-1);
        U = U(1:end-1);
        V = V(1:end-1);
    else
        n = length(P);
    end
    U(isnan(U)) = 0;
    V(isnan(V)) = 0;
    
    %Spectra settings
    nfft = floor(0.25*n);
    window = hanning(n,'periodic');
    min_f = 0.01; % mininum frequency, below which no correction is applied
    max_f = 1.5; % maximum frequency, above which no correction is applied
    
    g = 9.81;
    rho = 1.009629812357770e+03;
    for i = 1:length(p)
        g(i,:) = 9.780318*(1.0+(5.2788E-3+2.36E-5*x)*x)+1.092E-6*p(i,:);
        h(i,:) = ((((-1.82E-15*p(i,:)+2.279E-10)*p(i,:)-2.2512E-5)*p(i,:)+9.72659)*p(i,:))/g(i,:);
    end
    h = mean(h);
    ws = puv(P,U,V,h,zp,zuv,fs,nfft,rho,window,min_f,max_f);
    wvstats.vta.vp3.time(it,1) = time(idxx(it));
    wvstats.vta.vp3.wavev(it,1) = ws.ubr;
    wvstats.vta.vp3.azm(it,1) = ws.azr;
    it = it+1;
end
save(['c:\Users\bkn5\Projects\Mekong_W2015\DataAnalysis\Spectra\Paper1\' 'WaveStats_HTA_VTA'],'wvstats','-v7.3')

