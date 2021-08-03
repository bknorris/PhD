%This script calculates Wave Energy Flux (WEF) from pressure data

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
close all
file = 'd:\Projects\Mekong_W2015\Data\Vector\DPS2\V5109_130315.mat';
% file = 'd:\Projects\Mekong_W2015\Data\Aquadopp\DPS2\AD5116_15March2015.mat';
sdir = 'd:\Projects\Mekong_W2015\DataAnalysis\Paper1\WPF\';
load(file)
dname = 'VTA'; %deployment
name = '5116'; %instrument
zp = 0.105; %pressure sensor height above bottom
intv = 10; %step interval (seconds)
win = 300; %window (seconds)
start = datenum(2015,03,14,07,10,00);
stop = datenum(2015,03,14,09,25,00);
savedatfile = 0;

%%%Instrument type setup variables
if exist('aqdp','var')
    %for Aquadopps:
    ind = find(aqdp.datenum >= start & aqdp.datenum <=stop);
    pressure = aqdp.pressure(ind);
    temperature = aqdp.temperature(ind);
    time = aqdp.datenum(ind);
    fs = str2double(aqdp.metadata.samprate(1:2));
    x = (sin(aqdp.metadata.lat/57.29578))^2;
elseif exist('ADV','var')
    %for ADVs:
    ind = find(ADV.datetime >= start & ADV.datetime <=stop);
    pressure = ADV.Pres(ind);
    time = ADV.datetime(ind);
    ind2 = find(ADV.Sensor.Datetime >= start & ADV.Sensor.Datetime <=stop);
    temp = ADV.Sensor.Temp(ind2); %Vector records temperature at 1Hz
    temperature = spline(ADV.Sensor.Datetime(ind2),temp,ADV.datetime(ind));
    fs = ADV.Metadata.instmeta.samprate; %sampling frequency for ADVs
    x = (sin(ADV.Metadata.inst_lat/57.29578))^2; %for ADVs
end
%%%Initialize Global Variables
WPF = struct();
avt = fs*intv; %# samples in interval
nwin = fs*win;
swin = fs*30; %30 second averaging window (to smooth)
DOF = round((nwin/swin)*2);
disp(['50% Hamming windowed spectra with ' sprintf('%0.0f',DOF) ' degrees of freedom'])
lf = 1.2;
hf = 0.05; 
nsamp = length(pressure);
ind = [1 avt:avt:nsamp];
for ii = 1:length(ind)
    if abs(nsamp-ind(ii)) < nwin  %skip the last few indexes approaching the end of the t-s
        continue
    else
        idx = ind(ii):ind(ii)+nwin-1;
    end
    p = pressure(idx)+zp;
    t = temperature(idx);
    p = cmgbridge(p,100,100,1000);
    t = cmgbridge(t,100,100,1000);
    if ~isempty(find(isnan(p),1))
        p(isnan(p)) = nanmean(p);
    end
    if ~isempty(find(isnan(t),1))
        t(isnan(t)) = nanmean(t);
    end
    time2 = time(ind(ii));
    salt = repmat(20,length(t),1); %assume constant salinity of 20PSU
    rho = SeaDensity(salt,t,p);
    
    g = zeros(length(p),1);h = zeros(length(p),1);
    for j = 1:length(p)
        g(j,:) = 9.780318*(1.0+(5.2788E-3+2.36E-5*x)*x)+1.092E-6*p(j,:);
        h(j,:) = ((((-1.82E-15*p(j,:)+2.279E-10)*p(j,:)-2.2512E-5)*p(j,:)+9.72659)*p(j,:))/g(j,:);
    end
    rho = mean(rho);
    h = mean(h);
    g = mean(g);
    p = detrend(p);
    
    [Cpp,fq] = pwelch(p,hanning(swin),swin*0.5,swin,fs);
    lfc = find(fq >= hf,1,'first');hfc = find(fq <= lf,1,'last');
    df = fq(3)-fq(2);
    omega = 2*pi.*fq;
    k = qkhf(omega,h)./h;
    kh = k*h;
    kz = k*zp;
    attn = cosh(kz)./cosh(kh);
    attn(attn<0.2) = 0.2;
    Spp = Cpp./(attn.^2);           %surface elevation spectrum
    m0 = sum(Spp(lfc:hfc)*df);
    Hs = 4*sqrt(m0);                %sig. wave height
    Hrms = Hs/sqrt(2);
    
    c = sqrt(g.*tanh(kh))./sqrt(k);
    n = (1/2)*(1+((2*kh)./(sinh(2*kh))));
    Cg = n.*c;                      %group speed
    Cg = Cg(lfc:hfc);
    E = (1/8)*rho.*g.*Hrms.^2;      %wave energy
    F = sum(E.*Cg);                 %wave power
    if F < 0
        F = 0; %just in case
    end
    
    j = 1:length(Cg);
    WPF.Cg(ii,j) = Cg;
    WPF.freq(ii,:) = fq;
    WPF.E(ii,:) = E;
    WPF.F(ii,:) = F;
    WPF.h(ii,:) = h;
    WPF.Hs(ii,:) = Hs;
    WPF.time(ii,:) = time2;
end

WPF.Info.DepStart = [datestr(start,'dd-mm-yyyy HH:MM:SS') ' ICT'];
WPF.Info.DepStop = [datestr(stop,'dd-mm-yyyy HH:MM:SS') ' ICT'];
WPF.Info.avgint = [num2str(win) ' second averages'];
WPF.Info.step = [num2str(intv) ' second step between averages'];
WPF.Info.cmt = 'Each row represents a new interval';
WPF.Info.timecmt = 'Times listed are the user-specified start/stop times used to generate these data';

if savedatfile
    filename = [name '_' dname];
    sfname = ['Wpf' filename];
    save([sdir sfname],'WPF','-v7.3')
    disp(['WPF data file saved as ' sfname])
end






