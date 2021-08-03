%Cutoff frequencies for the wave/turbulence band using the method of
%Bricker and Monismith (2007)

clear
%%%Data Pre-Processing%%%
start = datenum('07-Mar-2015 14:10:00');
stop = datenum('07-Mar-2015 17:10:00');
inst = 'd:\Projects\Mekong_W2015\Data\Vector\FSS\VC01_070315.mat';
% start = datenum('08-Mar-2015 14:25:00');
% stop = datenum('08-Mar-2015 19:00:00');
% inst = 'd:\Projects\Mekong_W2015\Data\Vector\FSS\VC01_080315.mat';
% start = datenum('10-Mar-2015 14:45:00');
% stop = datenum('10-Mar-2015 16:30:00');
% inst = 'd:\Projects\Mekong_W2015\Data\Vector\FSS\V5109_100315.mat';
% start = datenum('13-Mar-2015 16:15:40');
% stop = datenum('13-Mar-2015 20:04:40');
% inst = 'd:\Projects\Mekong_W2015\Data\Aquadopp\DPS2\AD5116_15March2015.mat';
% start = datenum('14-Mar-2015 04:40:00');
% stop = datenum('14-Mar-2015 10:40:00');
% inst = 'd:\Projects\Mekong_W2015\Data\Aquadopp\DPS2\AD5116_15March2015.mat';
%Cutoffs
lf = 2;
hf = 0.05;
%Averaging interval
win = 600; %seconds (10 minutes)
step = 60; %seconds

%%%DATA PROCESSING%%%
disp(['Loading ' inst])
load(inst)
if exist('ADV','var')
    zp = ADV.Metadata.pressure_sensor_height/1000;
    zuv = ADV.Metadata.velocity_height/1000;
    time1 = ADV.datetime;
    vid = find(time1 >= start & time1 <= stop);
    time1 = time1(vid);
    U = ADV.U(vid);V = ADV.V(vid);W = ADV.W(vid);P = ADV.Pres(vid);
    x = (sin(ADV.Metadata.inst_lat/57.29578))^2;
    fs = ADV.Metadata.instmeta.samprate;
    avt = fs*step;
    nwin = fs*win;
    swin = fs*40; %30 second averaging window (to smooth)
    DOF = round((nwin/swin)*2);
    disp(['50% Hamming windowed spectra with ' sprintf('%0.0f',DOF) ' degrees of freedom'])
elseif exist('aqdp','var')
    zp = aqdp.metadata.HAB/1000;
    zuv = zp;
    time1 = aqdp.datenum;
    vid = find(time1 >= start & time1 <= stop);
    time1 = time1(vid);
    U = nanmean(aqdp.u(vid),2);V = nanmean(aqdp.v(vid),2);W = nanmean(aqdp.w(vid),2);P = nanmean(aqdp.pressure(vid),2);
    x = (sin(aqdp.metadata.lat/57.29578))^2;
    fs = str2double(regexprep(aqdp.metadata.samprate,' Hz',''));
    avt = fs*step;
    nwin = fs*win;
    swin = fs*40; %30 second averaging window (to smooth)
    DOF = round((nwin/swin)*2);
    disp(['50% Hamming windowed spectra with ' sprintf('%0.0f',DOF) ' degrees of freedom'])
end
U = cmgbridge(U,100,100,1000);V = cmgbridge(V,100,100,1000);
W = cmgbridge(W,100,100,1000);P = cmgbridge(P,100,100,1000);
U(isnan(U)) = 0;V(isnan(V)) = 0;W(isnan(W)) = 0;P(isnan(P)) = 0;

%%%Data Analysis
nsamp = length(vid);
ind = [1 avt:avt:nsamp];
disp('Starting analysis')
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
    time = time1(ind(ii));
    
    g = zeros(length(p),1);h = zeros(length(p),1);
    for j = 1:length(p)
        g(j,:) = 9.780318*(1.0+(5.2788E-3+2.36E-5*x)*x)+1.092E-6*p(j,:);
        h(j,:) = ((((-1.82E-15*p(j,:)+2.279E-10)*p(j,:)-2.2512E-5)*p(j,:)+9.72659)*p(j,:))/g(j,:);
    end
    h = mean(h);
    g = mean(g);
    p = detrend(p);
    w = detrend(w);
    nf = length(p);
    rho = 1024; %kg/m^3
    
    win = 0.54-0.46*cos(2*pi*[0:nf-1]'/(nf-1)); %hamming window
    pf=(fft(p));
    pf = pf(1:nf/2+1); %get fft(ps) corresponding to [0 nf)
    f = 0:fs/length(p):fs/2;
    omega = 2*pi.*f;
    k = qkhf(omega,h)./h;
    kh = k*h;
    kzp = k*zp;
    %Surface Pressure Scaling
    attn  = ones(length(f),1);
    attn(2:end) = (cosh(kzp(2:end)./cosh(kh(2:end))));
    attn(attn<0.2) = 0.2;
    pfa = pf./(attn.^2);
    sp = ifft(pfa,length(p),'symmetric'); %scaled surface pressure
    %find coherence between w and p
    [msc,f] = mscohere(sp,w,hamming(swin),swin*0.5,nwin,fs);
    [Cww,~] = pwelch(w,hamming(swin),swin*0.5,nwin,fs);
    chl = cohere_signif_level(DOF);
    band = find(f>=hf&f<=lf);
    if nnz(msc(band)>chl) == 0
        fc = NaN;
    else
        idx = find(msc(band)>chl,1,'last');
        fc = f(idx);
    end
    %save to structure
    wave.time(ii) = time;
    wave.fc(ii) = fc;
end
fprintf('Max cutoff frequency: %0.2f\n',nanmax(wave.fc))
fprintf('Mean cutoff frequency: %0.2f\n',nanmean(wave.fc))
disp('End of File')