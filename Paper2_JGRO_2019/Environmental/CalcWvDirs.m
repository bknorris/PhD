%Calculate wave directions for the HTA and VTA
clear
%times are from wavy sections of the respective deployments during the
%latter stage of the rising tide
start(1) = datenum(2015,03,07,13,36,32);stop(1) = datenum(2015,03,07,17,05,02);
start(2) = datenum(2015,03,08,14,14,22);stop(2) = datenum(2015,03,08,19,01,21);
start(3) = datenum(2015,03,10,14,45,06);stop(3) = datenum(2015,03,10,16,36,36);
start(4) = datenum(2015,03,07,13,36,32);stop(4) = datenum(2015,03,08,08,00,00);
start(5) = datenum(2015,03,08,13,00,00);stop(5) = datenum(2015,03,10,16,36,36);
start(6) = datenum(2015,03,14,07,00,00);stop(6) = datenum(2015,03,14,13,17,33);
start(7) = datenum(2015,03,14,07,00,00);stop(7) = datenum(2015,03,14,13,17,33);

data = struct(); %preallocate data structure
%%%Data Pre-Processing%%%
%process ADV data first; everything will be rotated to principal wave axis
inst{1} = 'D:\Projects\Mekong_W2015\Data\Vector\FSS\V5109_070315.mat';
inst{2} = 'D:\Projects\Mekong_W2015\Data\Vector\FSS\V5109_080315.mat';
inst{3} = 'D:\Projects\Mekong_W2015\Data\Vector\FSS\V5109_100315.mat';
inst{4} = 'D:\Projects\Mekong_W2015\Data\Vector\SW_Mudflat\V5108_030315.mat';
inst{5} = 'D:\Projects\Mekong_W2015\Data\Vector\SW_Mudflat\V5108_080315.mat';
inst{6} = 'D:\Projects\Mekong_W2015\Data\Aquadopp\DPS2\AD5116_15March2015.mat';
inst{7} = 'D:\Projects\Mekong_W2015\Data\Vector\DPS2\V5109_130315.mat';
name = {'HTA1';'HTA2';'HTA4';'SW1';'SW2';'VTA';'NE1'};
%%%Spectra settings
fn = fieldnames(data);
zp = [0.09 0.09 0.09 0.196 0.196 0.105 0.12];
zuv = [0.398 0.660 0.854 0.380 0.315 0.105 0.29];
lf = 1.2;
hf = 0.05;
win = 600; %seconds (9 minutes)
step = 300; %seconds
for i = 1:7
    disp(['Loading ' inst{i}])
    load(inst{i})
    if exist('ADV','var')
        time1 = ADV.datetime;
        vid = find(time1 >= start(i) & time1 <= stop(i));
        U = ADV.U(vid);V = ADV.V(vid);P = ADV.Pres(vid);
        x = (sin(ADV.Metadata.inst_lat/57.29578))^2;
        fs = 32;
        avt = fs*step;
        nwin = fs*win;
        swin = fs*40; %30 second averaging window (to smooth)
        DOF = round((nwin/swin)*2);
        disp(['50% Hamming windowed spectra with ' sprintf('%0.0f',DOF) ' degrees of freedom'])
    elseif exist('aqdp','var')
        time1 = aqdp.datenum;
        vid = find(time1 >= start(i) & time1 <= stop(i));
        U = nanmean(aqdp.u(vid),2);V = nanmean(aqdp.v(vid),2);P = nanmean(aqdp.pressure(vid),2);
        x = (sin(aqdp.metadata.lat/57.29578))^2;
        fs = 8;
        avt = fs*step;
        nwin = fs*win;
        swin = fs*50; %30 second averaging window (to smooth)
        DOF = round((nwin/swin)*2);
        disp(['50% Hamming windowed spectra with ' sprintf('%0.0f',DOF) ' degrees of freedom'])
    end
    U = cmgbridge(U,100,100,1000);V = cmgbridge(V,100,100,1000);
    P = cmgbridge(P,100,100,1000);
    U(isnan(U)) = 0;V(isnan(V)) = 0;P(isnan(P)) = 0;
    vname = lower(name{i});
    time1 = time1(vid);
    
%     save variables to structure
    data.(vname).p = P;
    data.(vname).u = U;
    data.(vname).v = V;
    data.(vname).time = time1;
    data.(vname).x = x;
    
    %%Data Analysis
    nsamp = length(vid);
    ind = [1 avt:avt:nsamp];
    for ii = 1:length(ind)
        if abs(nsamp-ind(ii)) < nwin  %skip the last few indexes approaching the end of the t-s
            continue
        else
            idx = ind(ii):ind(ii)+nwin-1;
        end
%         spectra of ADVs
        u = U(idx);
        v = V(idx);
        p = P(idx)+zp(i);
        time = time1(ind(ii));
        
%         calculate RMS velocities (Luhar et al. 2013)
        Ec = (1/nwin)*sum(u);Nc = (1/nwin)*sum(v);
        Ewrms = sqrt((1/nwin)*sum((u-Ec).^2));
        Nwrms = sqrt((1/nwin)*sum((v-Nc).^2));
        Uc = sqrt(Ec^2+Nc^2);Uwrms = sqrt(Ewrms^2+Nwrms^2);
        Uw = sqrt(2)*Uwrms;
        
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
%         sig wave height should be calculated with non-directional
%         spectra, i.e. from the pressure t-s
        [Cpp,F] = pwelch(p,hanning(swin),swin*0.5,swin,fs);
        [Cuu,~] = pwelch(u,hanning(swin),swin*0.5,swin,fs);
        [Cvv,~] = pwelch(v,hanning(swin),swin*0.5,swin,fs);
        [Cpu,~] = cpsd(p,u,hanning(swin),swin*0.5,swin,fs);
        [Cpv,~] = cpsd(p,v,hanning(swin),swin*0.5,swin,fs);
        [Cuv,~] = cpsd(u,v,hanning(swin),swin*0.5,swin,fs);
        Guv = Cuu+Cvv;
        
%         wave parameters
        lfc = find(F >= hf,1,'first');hfc = find(F <= lf,1,'last');
        df = F(3)-F(2);
        omega = 2*pi.*F;
        k = qkhf(omega,h)./h;
        kh = k*h;
        kz = k*zuv(i);
        attn = cosh(kz)./cosh(kh);
        attn(attn<0.2) = 0.2;
        Spp = Cpp./(attn.^2);           %surface elevation spectrum
        m0 = sum(Spp(lfc:hfc)*df);
        Hs = 4*sqrt(m0);                %sig. wave height
%         compute amplitude spectrum for direction/spreading
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
        
%         save variables to structure
        data.(vname).time2(ii) = time;
        data.(vname).Hs(ii) = Hs;
        data.(vname).Tr(ii) = Tr;
        data.(vname).F(ii,:) = F(lfc:hfc);
        data.(vname).Dir(ii,:) = Dir;
        data.(vname).Spread(ii,:) = Spread;
        data.(vname).omegar(ii) = omegar;
        data.(vname).Uc(ii) = Uc;
        data.(vname).Uwrms(ii) = Uwrms;
    end
    clear ADV aqdp
end
save('d:\Projects\Mekong_W2015\DataAnalysis\Paper2\WaveStats_Full\WVStats_HTA_VTAv2','data','-v7.3')
