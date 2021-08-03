%Calulates power spectra & cross spectra for Turbulence signals to
%calculate lag times between Vectrinos. Also calculates wave parameters:
%wave orbital velocity, wave excursion, Uinf_rms and Um

%calculate wave velocity along PCA - use VC01 instead of V5109 to calculate
%wave orb. velocity & wave excursion -  it's above the canopy!

clear
savefigdir = 'd:\Projects\Mekong_W2015\Figures\Paper1\WaveVelsTurbulence\';
plotfigures = 0; %[off on]
bin = 5;
data = struct(); %preallocate data structure
start = datenum(2015,03,07,15,30,00);stop = datenum(2015,03,07,16,00,00);

%%%Data Pre-Processing%%%
%process ADV data first; everything will be rotated to principal wave axis
vectors{1} = 'D:\Projects\Mekong_W2015\Data\Vector\FSS\VC01_070315.mat';
vectors{2} = 'D:\Projects\Mekong_W2015\Data\Vector\FSS\V5109_070315.mat';
vname = {'V1';'V5'};
for i = 1:2
    load(vectors{i})
    time1 = ADV.datetime;
    vid = find(time1 >= start & time1 <= stop);
    u = ADV.U(vid);v = ADV.V(vid);p = ADV.Pres(vid);
    u = cmgbridge(u,100,100,1000);v = cmgbridge(v,100,100,1000);
    p = cmgbridge(p,100,100,1000);
    x = (sin(ADV.Metadata.inst_lat/57.29578))^2;
    if i == 1
        cmp = pca(u,v,[],0);
%         heading = abs(360-cmp.mdir);
%         rot = (pi*heading)/180;
    end
%     T = [cos(rot) -sin(rot);...
%         sin(rot) cos(rot)];
%     vels = [u v];
%     V = vels*T';
%     u = V(:,1);v = V(:,2);
    data.(vname{i}).p = p;
    data.(vname{i}).u = u;
    data.(vname{i}).v = v;
    data.(vname{i}).time = time1(vid);
    data.(vname{i}).x = x;
    clear ADV
end

%load vectrino/turbulence files from HTA1
datdir = 'd:\Projects\Mekong_W2015\DataAnalysis\Paper1\Eps_Vel_Spectra\';
load([datdir 'HTA_07153000_30min.mat'])
load([datdir 'Turbulence\' 'HTA_07153000_30minTKE.mat']);
fn = fieldnames(dat);
for i = 2:length(fn) %only VP2 and VP3
    %average turbulence statistics together
    eps = (Stat.(fn{i}).beam1.E(bin,:)+Stat.(fn{i}).beam2.E(bin,:)...
        +Stat.(fn{i}).beam3.E(bin,:)+Stat.(fn{i}).beam4.E(bin,:))./4;
    eps = cmgbridge(eps,10,100,100);eps(isnan(eps)) = 0;
    time1 = Stat.(fn{i}).time;time2 = dat.(fn{i}).time;
    %interpolate turbulence (@ 1 Hz) up to 50 Hz
    eps = spline(time1,eps,time2);eps(eps < 0) = 0;
    %load x,y data from dat structure, rotate to principal wave axis
    x = dat.(fn{i}).x(:,bin);y = dat.(fn{i}).y(:,bin);
    heading = abs(360-cmp.mdir);
    rot = (pi*heading)/180;
    T = [cos(rot) -sin(rot);...
        sin(rot) cos(rot)];
    vels = [x y];
    V = vels*T';
    x = V(:,1);y = V(:,2);
    data.(fn{i}).eps = eps;
    data.(fn{i}).x = detrend(x);
    data.(fn{i}).y = detrend(y);
    data.(fn{i}).time = time2;
end
clear dat Stat

%interpolate Vector data to 50Hz
for i = 1:2
    time2 = data.vpro3.time;
    time1 = data.(vname{i}).time;
    u = spline(time1,data.(vname{i}).u,time2);
    v = spline(time1,data.(vname{i}).v,time2);
    p = spline(time1,data.(vname{i}).p,time2);
    time = spline(time1,data.(vname{i}).time,time2);
    data.(vname{i}).p = p;
    data.(vname{i}).u = detrend(u);
    data.(vname{i}).v = detrend(v);
    data.(vname{i}).time = time;
end
%%%END Data Pre-processing%%%

%%%Spectra%%%
fn = fieldnames(data);
zp = [0.985 0.24]; %Vector pressure HAB
zuv = [0.985 0.398];
fs = 50;
win = 300; %seconds (5 minutes)
step = 1; %seconds
avt = fs*step;
nwin = fs*win;
swin = fs*30; %30 second averaging window (to smooth)
DOF = round((nwin/swin)*2);
sl = cohere_signif_level(DOF,0.95);
lf = 1.2; %low-freq cutoff
hf = 0.05; %high-freq cutoff
wave = struct();
for i = 3
    nsamp = length(data.(fn{i}).time);
    ind = [1 avt:avt:nsamp];
    for ii = 1:length(ind)
        if abs(nsamp-ind(ii)) < nwin  %skip the last few indexes approaching the end of the t-s
                continue
        else
            idx = ind(ii):ind(ii)+nwin-1;
        end
        if i < 3
            %spectra of ADVs
            u = data.(vname{i}).u(idx);
            v = data.(vname{i}).v(idx);
            p = data.(vname{i}).p(idx)+zp(i);
            x = data.(vname{i}).x;
            time = data.(vname{i}).time(ind(ii));
            
            %calculate RMS velocities (Luhar et al. 2013)
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
            %sig wave height should be calculated with non-directional
            %spectra, i.e. from the pressure t-s
            [Cpp,F] = pwelch(p,hanning(swin),swin*0.5,swin,fs);
            [Cuu,~] = pwelch(u,hanning(swin),swin*0.5,swin,fs);
            [Cvv,~] = pwelch(v,hanning(swin),swin*0.5,swin,fs);
            Cuv = Cuu+Cvv;
            
            %wave parameters
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
%             t = 0:(1/fs):win;
%             f = F(lfc:hfc);
%             for r = 1:length(f)
%                 w = 2*pi/f(r);
%                 k = qkhf(w,h)./h;
%                 kh = k*h;
%                 kz = k*zuv(i);
%                 Uinfmax = a*w*(cosh(kz)/sinh(kh));
%                 for q = 1:length(t)
%                     Uinfw(q) = Uinfmax.*cos(w.*t(q));
%                 end
%                 U2(r) = sum(Uinfw.^2);
%             end
%             pos= find(u>0); %onshore vels only
%             Uinfw = u(pos);
%             
%             T = max(1./F(lfc:hfc));
%             Uinfrms = sqrt((1/T)*trapz(Uinfw.^2)*(1/fs));
                
            %parameters of interest
            omegar = sum(omega(lfc:hfc).*Cuv(lfc:hfc)*df)./...
                sum(Cuv(lfc:hfc)*df);
            Tr = 2*pi/omegar;
%             T = 1./F(lfc:hfc);
%             Uw = a*omega.*(cosh(kz)./sinh(kh));
%             Uinfw = Uw(lfc:hfc).*cos(omega(lfc:hfc));
%             Uinfrms = sqrt((1./T).*trapz(Uinfw.^2).*(1/50));
            if i == 1
%                 Uinf_max = max(Uw);
                %save variables to structure
                wave.(fn{i}).time(ii) = time;
                wave.(fn{i}).Hs(ii) = Hs;
                wave.(fn{i}).Tr(ii) = Tr;
                wave.(fn{i}).omegar(ii) = omegar;
                wave.(fn{i}).Uc(ii) = Uc;
                wave.(fn{i}).Uwrms(ii) = Uwrms;
            elseif i == 2
                Um = Uc+Uw;
                %save variables to structure
                wave.(fn{i}).time(ii) = time;
                wave.(fn{i}).Tr(ii) = Tr;
                wave.(fn{i}).omegar(ii) = omegar;
                wave.(fn{i}).Um(ii) = Um;
                wave.(fn{i}).Uc(ii) = Uc;
                wave.(fn{i}).Uwrms(ii) = Uwrms;
            end
        elseif i == 3
            %cospectra of turbulence
            eps1 = data.(fn{i}).eps(idx);
            eps2 = data.(fn{i+1}).eps(idx);
            time = data.(fn{i}).time(ind(ii));
            
            [Cuv,f] = cpsd(eps2,eps1,hanning(swin),swin*0.5,swin,fs);
            [MSC,~] = mscohere(eps2,eps1,hanning(swin),swin*0.5,swin,fs);
            lfc = find(f < 0.05); %low-freq cutoff
            Cuv(lfc) = NaN;MSC(lfc) = NaN;
            [~,id] = max(Cuv);
            phase = -angle(Cuv);theta = rad2deg(phase(id));
            F = f(id);
            deltaT = theta/(360*F);
            coh = MSC(id);
            if coh < sl
                %                 disp(['Cospectra is not significant at timestep ' datestr(time)])
                deltaT = NaN;
                coh = NaN;
            end
            wave.vps.time(ii) = time;
            wave.vps.deltaT(ii) = deltaT;
            wave.vps.coh(ii) = coh;
            wave.vps.siglvl(ii) = sl;
        end
    end
end




if plotfigures
    f1 = figure(1);
    set(f1,'PaperOrientation','portrait',...
        'position',[400 200   1000   500]);
    set(gcf,'color','w','paperpositionmode','auto')
    sp(1) = subplot(211);
    plot(time2,x2,'-m','linewidth',1.5),hold on
    plot(time2,x1,'-k','linewidth',1.5)
    set(gca,'xticklabel',[],'xlim',[time2(1) time2(end)])
    title('Cross-Shore velocity, HTA1 VP2 & VP3')
    sp(2) = subplot(212);
    p(2) = plot(time2,vp3,'-m','linewidth',1.5);hold on
    p(1) = plot(time2,vp2,'-k','linewidth',1.5);
    set(gca,'ylim',[0 0.03],...
        'xlim',[time2(1) time2(end)])
    title('Turbulence, HTA1 VP2 & VP3')
    legend(p,{'VP2';'VP3'});
    datetickzoom('x','HH:MM:SS','keepticks','keeplimits')
    ylabel('\epsilon (Wkg^-^1)')
    xlabel(['Time on ' datestr(time2(1),'dd-mm-yyyy')])
    linkaxes(sp,'x')
    prettyfigures('font','arial','text',11,'labels',13,'box',1,'grid',1)
end


