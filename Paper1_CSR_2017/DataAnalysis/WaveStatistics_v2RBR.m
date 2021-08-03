%Wave Statistics for RBRs (Duet only): H, Hs, T
start = datenum(2015,03,14,07,00,00);
stop = datenum(2015,03,14,13,19,23);

%%%Data Pre-Processing%%%
inst = 'D:\Mekong_W2015\Data\RBR\Duet\DPS2\Duet_140315.mat';
%Cutoffs
lf = 1.2;
hf = 0.05;
%Averaging interval
win = 60; %seconds (1 minutes)
step = 30; %seconds

%%%DATA PROCESSING%%%
tic
disp(['Loading ' inst])
load(inst)
zp = str2double(regexprep(RBR.Metadata.hab,' mm',''))/1000;
time1 = RBR.Datetime;
vid = find(time1 >= start & time1 <= stop);
time1 = time1(vid);
P = nanmean(RBR.SeaPres(vid),2);
x = (sin(RBR.Metadata.latitude/57.29578))^2;
fs = str2double(regexprep(RBR.Metadata.samprate,'Hz',''));
avt = fs*step;
nwin = fs*win;
swin = fs*5; %10 second averaging window (to smooth)
DOF = round((nwin/swin)*2);
disp(['50% Hamming windowed spectra with ' sprintf('%0.0f',DOF) ' degrees of freedom'])
P = cmgbridge(P,100,100,1000);
P(isnan(P)) = 0;

%save variables to structure
wave.p = P;
wave.time = time1;
wave.x = x;

%%%Data Analysis
nsamp = length(vid);
ind = [1 avt:avt:nsamp];
for ii = 1:length(ind)
    if abs(nsamp-ind(ii)) < nwin  %skip the last few indexes approaching the end of the t-s
        continue
    else
        idx = ind(ii):ind(ii)+nwin-1;
    end
    %spectra of ADVs
    p = P(idx)+zp;
    time = time1(ind(ii));
    
    
    g = zeros(length(p),1);h = zeros(length(p),1);
    for j = 1:length(p)
        g(j,:) = 9.780318*(1.0+(5.2788E-3+2.36E-5*x)*x)+1.092E-6*p(j,:);
        h(j,:) = ((((-1.82E-15*p(j,:)+2.279E-10)*p(j,:)-2.2512E-5)*p(j,:)+9.72659)*p(j,:))/g(j,:);
    end
    h = mean(h);
    g = mean(g);
    p = detrend(p);
    %sig wave height should be calculated with non-directional
    %spectra, i.e. from the pressure t-s
    [Cpp,F] = pwelch(p,hanning(swin),swin*0.5,swin,fs);
    
    %wave parameters
    lfc = find(F >= hf,1,'first');hfc = find(F <= lf,1,'last');
    df = F(3)-F(2);
    omega = 2*pi.*F;
    k = qkhf(omega,h)./h;
    kh = k*h;
    kz = k*zp;
    attn = cosh(kz)./cosh(kh);
    attn(attn<0.2) = 0.2;
    Spp = Cpp./(attn.^2);           %surface elevation spectrum
    m0 = sum(Spp(lfc:hfc)*df);
    Hs = 4*sqrt(m0);                %sig. wave height
    [pk,ll] = findpeaks(Cpp);
    [~,id] = max(pk);
    T = 1./F(ll(id));
    
    %save variables to structure
    wave.time2(ii) = time;
    wave.h(ii) = h;
    wave.Hs(ii) = Hs;
    wave.T(ii) = T;
    wave.F(ii,:) = F(lfc:hfc);
end
disp(['Analysis completed in: ' num2str(toc/60) ' minutes'])

%%%Save Data%%%
sid = strfind(inst,'\Data\');
path1 = [inst(1:sid) 'DataAnalysis\'];
name = regexp(inst,'\.*\w\','split');name = regexprep(name{end},'.mat','');
save([path1 name 'VTA2wvs'],'wave','-v7.3')

