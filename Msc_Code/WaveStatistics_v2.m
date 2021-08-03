%Wave Statistics: Hs, T, Dir, Spread, orbital v, Uc, Uwrms for any Aquadopp
%or ADV

clear
start = datenum('07-Mar-2015 14:00:00');
stop = datenum('07-Mar-2015 17:00:00');

%%%Data Pre-Processing%%%
inst = 'd:\Projects\Mekong_W2015\Data\Vector\FSS\V5109_070315.mat';
%Cutoffs
lf = 1.2;
hf = 0.05;
%Averaging interval
win = 180; %seconds (5 minutes)
step = 30; %seconds

%%%DATA PROCESSING%%%
tic
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
    swin = fs*10; %30 second averaging window (to smooth)
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
    swin = fs*10; %10 second averaging window (to smooth)
    DOF = round((nwin/swin)*2);
    disp(['50% Hamming windowed spectra with ' sprintf('%0.0f',DOF) ' degrees of freedom'])
end
U = cmgbridge(U,100,100,1000);V = cmgbridge(V,100,100,1000);
W = cmgbridge(W,100,100,1000);P = cmgbridge(P,100,100,1000);
U(isnan(U)) = 0;V(isnan(V)) = 0;W(isnan(W)) = 0;P(isnan(P)) = 0;

%save variables to structure
wave.p = P;
wave.u = U;
wave.v = V;
wave.w = W;
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
    u = U(idx);
    v = V(idx);
    w = W(idx);
    p = P(idx)+zp;
    time = time1(ind(ii));
    
    %calculate RMS velocities (Luhar et al. 2013)
    Ec = (1/nwin)*sum(u);Nc = (1/nwin)*sum(v);
    Ewrms = sqrt((1/nwin)*sum((u-Ec).^2));
    Nwrms = sqrt((1/nwin)*sum((v-Nc).^2));
    Wc = (1/nwin)*sum(w);Wwrms = sqrt((1/nwin)*sum((w-Wc).^2));
    Uc = sqrt(Ec^2+Nc^2);Uwrms = sqrt(Ewrms^2+Nwrms^2);
%     Uw = sqrt(2)*Uwrms;
%     Ww = sqrt(2)*Wwrms;
    
    g = zeros(length(p),1);h = zeros(length(p),1);
    for j = 1:length(p)
        g(j,:) = 9.780318*(1.0+(5.2788E-3+2.36E-5*x)*x)+1.092E-6*p(j,:);
        h(j,:) = ((((-1.82E-15*p(j,:)+2.279E-10)*p(j,:)-2.2512E-5)*p(j,:)+9.72659)*p(j,:))/g(j,:);
    end
    h = mean(h);
    g = mean(g);
    p = detrend(p);
    u = detrend(u);
    v = detrend(v);
    w = detrend(w);
    %sig wave height should be calculated with non-directional
    %spectra, i.e. from the pressure t-s
    [Cpp,F] = pwelch(p,hanning(swin),swin*0.5,swin,fs);
    [Cuu,~] = pwelch(u,hanning(swin),swin*0.5,swin,fs);
    [Cvv,~] = pwelch(v,hanning(swin),swin*0.5,swin,fs);
    [Cww,~] = pwelch(w,hanning(swin),swin*0.5,swin,fs);
    [Cpu,~] = cpsd(p,u,hanning(swin),swin*0.5,swin,fs);
    [Cpv,~] = cpsd(p,v,hanning(swin),swin*0.5,swin,fs);
    [Cuv,~] = cpsd(u,v,hanning(swin),swin*0.5,swin,fs);
    Guv = Cuu+Cvv;
    
    %wave parameters
    lfc = find(F >= hf,1,'first');hfc = find(F <= lf,1,'last');
    df = F(3)-F(2);
    omega = 2*pi.*F;
    k = qkhf(omega,h)./h;
    kh = k*h;
    kz = k*zuv;
    attn = cosh(kz)./cosh(kh);
    attn(attn<0.2) = 0.2;
    Spp = Cpp./(attn.^2);           %surface elevation spectrum
    m0 = sum(Spp(lfc:hfc)*df);
    Hs = 4*sqrt(m0);                %sig. wave height
    %compute amplitude spectrum for direction/spreading
    Cpu = Cpu.*conj(Cpu);Cpv = Cpv.*conj(Cpv);
    Cuu = Cuu.*conj(Cuu);Cvv = Cvv.*conj(Cvv);
    Cuv = Cuv.*conj(Cuv);
    Dir=57.296*atan2(Cpu(lfc:hfc),Cpv(lfc:hfc));
    Dir=mod(Dir+180,360)-180;
    R2=((Cuu(lfc:hfc)-Cvv(lfc:hfc)).^2+4*Cuv(lfc:hfc).^2).^.5./(Cuu(lfc:hfc)+Cvv(lfc:hfc));
    Spread = 57.296*((1-R2)/2).^.5; %wave spreading
    omegar = sum(omega(lfc:hfc).*Guv(lfc:hfc)*df)./...
        sum(Guv(lfc:hfc)*df);
    Tr = 2*pi/omegar;
    ubr = sqrt(2*sum(Guv(lfc:hfc)*df));
    %compute vertical velocity spectra (e.g. Gerbi et al. 2008)
%     Sww = Cpp.*((k.^2)./((1.025^2)*(omega.^2))).*(tanh(kh).*tanh(kz));
%     if nnz(u == 0)/length(u) > 0.85 %if 85% of burst is 0 ignore
%         kc = NaN;
%     else
%         [F0,~] = intersections(F,0.3.*Cww,F,Sww);
%         if isempty(F0)
%             kc = NaN;
%         else
%             kc = F0(1);
%         end
%     end
    
    %save variables to structure
    wave.time2(ii) = time;
    wave.h(ii) = h;
    wave.Hs(ii) = Hs;
    wave.Tr(ii) = Tr;
    wave.F(ii,:) = F(lfc:hfc);
    wave.Dir(ii,:) = Dir;
    wave.Spread(ii,:) = Spread;
    wave.ubr(ii) = omegar;
    wave.omegar(ii) = omegar;
    wave.Uc(ii) = Uc;
    wave.Uwrms(ii) = Uwrms;
    wave.Wc(ii) = Wc;
    wave.Wwrms(ii) = Wwrms;
%     wave.kc(ii) = kc; %wave cutoff frequency
end
disp(['Analysis completed in: ' num2str(toc/60) ' minutes'])
%%%Save Data%%%
% sid = strfind(inst,'\Data\');
% path1 = [inst(1:sid) 'DataAnalysis\Paper2\WaveStats_Full\'];
% name = regexp(inst,'\.*\w\','split');name = regexprep(name{end},'.mat','');
% save([path1 name 'wvs'],'wave','-v7.3')

