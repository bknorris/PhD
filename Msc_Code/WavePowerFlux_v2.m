%This script calculates Wave Power Flux (WPF) from pressure data

%%Preliminary Checklist:
%Terms (units):
%k = wave number
%z = instrument hab (m) [nominal height of pressure sensor]
%h = water depth (m)
%Hs = sig. wave height (m)
%omega = wave angular velocity (rad*m/s)
%g = gravity (m/s^2)
%rho = water density (kg/m^3)
%Cg = wave group velocity (m/s)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear %a good idea
load('HR3_9March2015_f_pad.mat')
dname = 'FSS32'; %deployment
name = 'ADHR3'; %instrument
ifaqdp = 1; %turn on if aquadopp
ifvector = 0; %turn on if vector
plotfigs = 0;
z = 0.075; %pressure sensor height above bottom
intv = 1; %averaging interval (minutes)
start = datenum(2015,03,08,15,45,06);
stop = datenum(2015,03,08,16,42,38);
savedatfile = 0;
% plot(aqdp.datenum,aqdp.pressure)
% datetickzoom('x','keepticks','keeplimits')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Instrument type setup variables
if ifaqdp
    %for Aquadopps:
    ind = find(aqdp.datenum >= start & aqdp.datenum <=stop);
    pressure = aqdp.pressure(ind);
    temperature = aqdp.temperature(ind);
    fs = str2double(aqdp.metadata.samprate(1:2));
    x = (sin(aqdp.metadata.lat/57.29578))^2;
elseif ifvector
    %for ADVs:
    ind = find(ADV.datetime >= start & ADV.datetime <=stop);
    pressure = ADV.Pres(ind);
    ind2 = find(ADV.Sensor.Datetime >= start & ADV.Sensor.Datetime <=stop);
    temp = ADV.Sensor.Temp(ind2); %Vector records temperature at 1Hz
    temperature = spline(ADV.Sensor.Datetime(ind2),temp,ADV.datetime(ind));
    fs = ADV.Metadata.instmeta.samprate; %sampling frequency for ADVs
    x = (sin(ADV.Metadata.inst_lat/57.29578))^2; %for ADVs
end
%Initialize Global Variables %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
WPF = struct();
avt = fs*intv*60; %# samples in interval
idx = avt:avt:length(pressure);
idxx = [1 idx];
it = 1;
%Iterate through time intervals %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
while it < length(idxx)
    disp(['Iteration number ' num2str(it)])
    ind = idxx(it):idxx(it+1);
    p = pressure(ind);
    p = p+z; %adjust for height of the pressure sensor
    t = temperature(ind);
    if rem(length(p), 2) == 1
        p = p(1:end-1); %force inputs to be even
        t = t(1:end-1);
    end
    %Fill gaps in temperature/pressure signal %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if cmgidgaps(p) > 0
        disp(['Found ' num2str(cmgidgaps(p)) ' gaps in pressure time-series'])
        disp(['Found ' num2str(cmgidgaps(t)) ' gaps in temperature time-series'])
        nlin = 1E2;
        maxgaps = 1E5;
        p = cmgbridge(p,nlin,maxgaps,maxgaps);
        t = cmgbridge(t,nlin,maxgaps,maxgaps);
        if cmgidgaps(p) > 0 || cmgidgaps(t) > 0
            p(isnan(p)) = nanmean(p);
            t(isnan(t)) = nanmean(t);
            disp(['Number of gaps in pressure remaining: ' num2str(cmgidgaps(p))])
            disp(['Number of gaps in temperature remaining: ' num2str(cmgidgaps(t))])
        else
            disp('Gaps filled')
        end
    end
    salt = repmat(20,length(t),1); %assume constant salinity of 20PSU
    rho = SeaDensity(salt,t,p);
    
    %Calculate Depth (h) from the pressure signal %%%%%%%%%%%%%%%%%%%%%%%%%
    %designate start/stop time
    %From: UNESCO (1983): Algorithms for computation of fundamental properties
    %of seawater. UNESCO technical papers in marine science 44:1-55.
    g = zeros(length(p),1);h = zeros(length(p),1);
    for i = 1:length(p)
        g(i,:) = 9.780318*(1.0+(5.2788E-3+2.36E-5*x)*x)+1.092E-6*p(i,:);
        h(i,:) = ((((-1.82E-15*p(i,:)+2.279E-10)*p(i,:)-2.2512E-5)*p(i,:)+9.72659)*p(i,:))/g(i,:);
    end
    h = mean(h);
    g = mean(g);
    
    %Calculate PSD(timeseries) %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ps = detrend(p);
    n = length(ps);
    rho = mean(rho);     
    nf = fs/2; %nyquist criterium
    nfft = n/2;
    min_f = 0.05; % mininum frequency, below which no correction is applied
    max_f = 0.33; % maximum frequency, above which no correction is applied
    [psd,freq] = pwelch(ps,hanning(n),n/2,n,fs,'psd');
    pf = psd(2:end,1);
    freq(1) = [];
    
    %Calculate wavenumber (k) and attenuation%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    disp('Computing surface pressure attenuation adjustment')
    k = zeros(length(freq),1);
    attn = zeros(length(freq),1); %attenuation
    T = 1./freq; %period = 1/freq
    omega = 2*pi./T';
    for i = 1:length(freq)
        k(i) = om2k(omega(i),h);
        %From Raubenheimer, B et al. 1998 "Estimating Wave Heights from
        %Pressure Measured in Sand Bed"
        attn(i) = (cosh(k(i)*(z)))/(cosh(k(i)*h));
    end
    attn(freq<min_f | freq>max_f) = 1;	% only for selected range
    attn(attn<0.2) = 0.2; % correction factor never higher than 5sec
    fmax=max(find(freq<=max_f));
	flim=fmax:min(fmax+fix(length(k)/10),length(k));
	attn(flim)=(length(flim):-1:1)*(attn(fmax)-1)/length(flim)+1; 

    %Calculate attenuated surface pressure %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    disp('Applying adjustment to timeseries')
    pfa = pf./(attn.^2); %adjusted fourier transfm.

    %Calculate Wave Parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    wstats = wave_stats(pfa,freq,max_f,ps,h,fs,z,n/2);
    %check Hsig
    df = freq(1);
    m0 = freq.^0.*pfa*df;
    Hs = 4*sqrt(m0);
    Hrms = Hs./sqrt(2);
    
    %Calculate wave group velocity %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %dfn. for Cg from Dalrymple et al. 1984 "Wave Diffraction due to Areas"
    kh = k.*h;
    c = sqrt(g.*tanh(kh))./sqrt(k);
    n = (1/2)*(1+((2*kh)./(sinh(2*kh))));
    Cg = n.*c;
    cgnan = find(isnan(Cg));
    Cg(cgnan) = Cg(cgnan+1);
    
    %Calculate energy (E) %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    E = (1/8)*rho.*g.*Hrms.^2;
    
    %Energy flux (F) %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    F = sum(E.*Cg);
    
    %Copy data to structure %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    j = 1:length(Cg);
    WPF.Cg(it,j) = Cg;
    WPF.E(it,:) = E;
    WPF.F(it,:) = F;
    WPF.h(it,:) = h;
    WPF.rho(it,:) = rho;
    WPF.freq(it,j) = freq;
    it = it+1;
end
disp('Finished')
figure
time = 1:length(WPF.F);
pf = polyfit(time',WPF.F,2);
pv = polyval(pf,time);
plot(time,WPF.F,'--or')
hold on
plot(time,pv,'Color','k','LineWidth',1.5)
xlabel('Time Elapsed (min)')
ylabel('Wave Power Flux (m^3/s^3)')
title('WPF Generation')
t = 10:10:length(WPF.F);t = [1 t];
for i = 1:length(t)-1
    wpf_linfit(i,:) = mean(pv(t(i):t(i+1)));
end

WPF.Ffit = wpf_linfit;
WPF.Info.DepStart = [datestr(start,'dd-mm-yyyy HH:MM:SS') ' ICT'];
WPF.Info.DepStop = [datestr(stop,'dd-mm-yyyy HH:MM:SS') ' ICT'];
WPF.Info.avgint = [num2str(intv) ' minute averages'];
WPF.Info.cmt = 'Each row represents a new interval';
WPF.Info.timecmt = 'Times listed are the user-specified start/stop times used to generate these data';

if savedatfile
    %save file
    filename = [name '_' dname];
    sfname = ['Wpf' filename];
    save(sfname,'WPF','-v7.3')
    disp(['Drag force data file saved as ' sfname])
end

