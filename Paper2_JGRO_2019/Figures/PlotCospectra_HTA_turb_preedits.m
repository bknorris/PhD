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
%     if i == 1
%         cmp = pca(u,v,[],0);
%         %         heading = abs(360-cmp.mdir);
%         %         rot = (pi*heading)/180;
%     end
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
for i = 1:length(fn) %only VP2 and VP3
    %average turbulence statistics together
    eps = mean((Stat.(fn{i}).beam1.E(1:bin,:)+Stat.(fn{i}).beam2.E(1:bin,:)...
        +Stat.(fn{i}).beam3.E(1:bin,:)+Stat.(fn{i}).beam4.E(1:bin,:))./4);
    eps = cmgbridge(eps,10,100,100);eps(isnan(eps)) = 0;
    time1 = Stat.(fn{i}).time;time2 = dat.(fn{i}).time;
    %interpolate turbulence (@ 1 Hz) up to 50 Hz
    eps = spline(time1,eps,time2);eps(eps < 0) = 0;
    %load x,y data from dat structure, rotate to principal wave axis
    x = mean(dat.(fn{i}).x(:,1:bin),2);y = mean(dat.(fn{i}).y(:,1:bin),2);
    %     heading = abs(360-cmp.mdir);
    %     rot = (pi*heading)/180;
    %     T = [cos(rot) -sin(rot);...
    %         sin(rot) cos(rot)];
    %     vels = [x y];
    %     V = vels*T';
    %     x = V(:,1);y = V(:,2);
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
win = 180; %seconds (5 minutes)
step = 10; %seconds
avt = fs*step;
nwin = fs*win;
swin = fs*10; %30 second averaging window (to smooth)
DOF = round((nwin/swin)*2);
disp(['50% Hamming windowed spectra with ' sprintf('%0.0f',DOF) ' degrees of freedom'])
lf = 1.2;
hf = 0.05; 
wave = struct();
for i = 1:5
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
            p = detrend(p);
            u = detrend(u);
            v = detrend(v);
            %sig wave height should be calculated with non-directional
            %spectra, i.e. from the pressure t-s
            [Cpp,F] = pwelch(p,hanning(swin),swin*0.5,swin,fs);
            [Cuu,~] = pwelch(u,hanning(swin),swin*0.5,swin,fs);
            [Cvv,~] = pwelch(v,hanning(swin),swin*0.5,swin,fs);
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
            kz = k*zuv(i);
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
            if i == 1
                %save variables to structure
                wave.(fn{i}).time(ii) = time;
                wave.(fn{i}).Hs(ii) = Hs;
                wave.(fn{i}).Tr(ii) = Tr;
                wave.(fn{i}).F(ii,:) = F(lfc:hfc);
                wave.(fn{i}).Dir(ii,:) = Dir;
                wave.(fn{i}).Spread(ii,:) = Spread;
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
        elseif i > 2
            %cross-correlation of turbulence
            %xcor is better than cpsd for non-sinusoid timeseries
            eps = data.(fn{i}).eps(idx);
            v = data.(fn{i}).x(idx);
            u = data.(fn{i}).y(idx);
            umag = sin(v).*sqrt(u.^2+v.^2); 
            time = data.(fn{i}).time(ind(ii));
            %split t-s into shorter windows (per wavelength) to increase
            %resolution
            [pf,f] = pwelch(umag,hanning(swin),swin*0.5,swin,fs);
            [~,id] = max(pf);
            T = 1/f(id);
            wvsamp = round(T*fs);
            idx2 = [1 wvsamp:wvsamp:length(umag)];
            timeDiff = zeros(length(idx2),1);
            for j = 1:length(idx2)-1
                [acor,lag] = xcorr(umag(idx2(j):idx2(j+1)),eps(idx2(j):idx2(j+1)));
                [~,I] = max(abs(acor));
                lagDiff = lag(I);
                timeDiff(j) = lagDiff/fs;
            end
            wave.(fn{i}).time(ii) = time;
            wave.(fn{i}).deltaT(ii) = var(timeDiff);
        end
    end
end
Ainfrms = wave.V1.Uwrms./wave.V1.omegar;

if plotfigures
    start = datenum(2015,03,07,15,45,00);
    stop = datenum(2015,03,07,15,46,00);
    f1 = figure(1);
    set(f1,'PaperOrientation','portrait',...
        'position',[400 200   1000   800]);
    set(gcf,'color','w','paperpositionmode','auto')
    fn = fieldnames(data); 
    time = data.vpro1.time;
    symb = {'o';'d';'^'};
    line = {'-';'--';'-.'};
    c = flipud([0.1 0.1 0.1;0.5 0.5 0.5;0.7 0.7 0.7]);
    ct = 1;
    sp(1) = subplot(411);
    for i = 3:5
        u = data.(fn{i}).y;
        plot(time,u,line{ct},...
            'color',c(ct,:),'linewidth',1.5),hold on
%         markx = time(1:10:end);
%         marky = u(1:10:end);
%         plot(markx,marky,symb{ct},...
%         'markersize',4,...
%         'markerfacecolor',c(ct,:),...
%         'markeredgecolor','k',...
%         'linewidth',1)
    ct = ct+1;
    end
    set(gca,'xlim',[start stop],'xticklabel',[])
    title('u velocities')
    sp(2) = subplot(412);
    ct = 1;
    for i = 3:5
        v = data.(fn{i}).x;
        plot(time,v,line{ct},...
            'color',c(ct,:),'linewidth',1.5),hold on
%         markx = time(1:10:end);
%         marky = v(1:10:end);
%         plot(markx,marky,symb{ct},...
%         'markersize',4,...
%         'markerfacecolor',c(ct,:),...
%         'markeredgecolor','k',...
%         'linewidth',1)
    ct = ct+1;
    end
    set(gca,'xlim',[start stop],'xticklabel',[])
    title('v velocities')
    sp(3) = subplot(413);
    ct = 1;
    for i = 3:5
        eps = data.(fn{i}).eps;
        plot(time,eps,line{ct},...
            'color',c(ct,:),'linewidth',1.5),hold on
%         markx = time(1:10:end);
%         marky = eps(1:10:end);
%         plot(markx,marky,symb{ct},...
%         'markersize',4,...
%         'markerfacecolor',c(ct,:),...
%         'markeredgecolor','k',...
%         'linewidth',1)
    ct = ct+1;
    end
    set(gca,'xlim',[start stop],'xticklabel',[])
    title('v velocities')
    sp(4) = subplot(414);
    ct = 1;
    time = wave.vpro1.time;
    pp = zeros(1,3);
    for i = 3:5
        dT = wave.(fn{i}).deltaT;
        pp(ct) = plot(time,dT,line{ct},...
            'color',c(ct,:),'linewidth',1.5);hold on
%         markx = time(1:10:end);
%         marky = dT(1:10:end);
%         plot(markx,marky,symb{ct},...
%         'markersize',4,...
%         'markerfacecolor',c(ct,:),...
%         'markeredgecolor','k',...
%         'linewidth',1)
    ct = ct+1;
    end
    set(gca,'xlim',[start stop],'ylim',[0 2.5])
    leg = legend(pp,{'x = -10 cm';'x = 10 cm';'x = 20 cm'});
    title('delta t')
    datetickzoom('x','HH:MM:SS','keepticks','keeplimits')
    linkaxes(sp,'x')


%     f2 = figure(2);
%     set(f2,'PaperOrientation','portrait',...
%         'position',[400 200   1000   500]);
%     ct = 1;
%    for i = 3:5 %VPs only
%     sp(1) = subplot(121);
%     plot(wave.(fn{i}).deltaT,wave.V1.Tr,symb{ct},...
%         'markerfacecolor',c(ct,:),...
%         'markeredgecolor','k','linewidth',1.5), hold on
%     ylabel('Peak Period (s)')
%     xlabel('\epsilon time lag (s)')
%     sp(2) = subplot(122);
%     plot(wave.(fn{i}).deltaT,wave.V1.Hs,symb{ct},...
%         'markerfacecolor',c(ct,:),...
%         'markeredgecolor','k','linewidth',1.5), hold on
%     ylabel('H_s (m)')
%     xlabel('\epsilon time lag (s)')
%     ct = ct+1;
%    end
    prettyfigures('font','arial','text',11,'labels',13,'tickdir','in','box',1,'grid',0)
    
        handles = findobj('type','figure');h = sort(handles);
    %     export_fig(h(1),[savefigdir 'Hs_Tr_epsDeltaT_ts'],'-png')
        export_fig(h(1),[savefigdir 'Hs_Tr_waveDeltaT'],'-png')
    
    %     f1 = figure(1);
    %     set(f1,'PaperOrientation','portrait',...
    %         'position',[400 200   1000   500]);
    %     set(gcf,'color','w','paperpositionmode','auto')
    %     sp(1) = subplot(211);
    %     plot(time2,x2,'-m','linewidth',1.5),hold on
    %     plot(time2,x1,'-k','linewidth',1.5)
    %     set(gca,'xticklabel',[],'xlim',[time2(1) time2(end)])
    %     title('Cross-Shore velocity, HTA1 VP2 & VP3')
    %     sp(2) = subplot(212);
    %     p(2) = plot(time2,vp3,'-m','linewidth',1.5);hold on
    %     p(1) = plot(time2,vp2,'-k','linewidth',1.5);
    %     set(gca,'ylim',[0 0.03],...
    %         'xlim',[time2(1) time2(end)])
    %     title('Turbulence, HTA1 VP2 & VP3')
    %     legend(p,{'VP2';'VP3'});
    %     datetickzoom('x','HH:MM:SS','keepticks','keeplimits')
    %     ylabel('\epsilon (Wkg^-^1)')
    %     xlabel(['Time on ' datestr(time2(1),'dd-mm-yyyy')])
    %     linkaxes(sp,'x')
    %     prettyfigures('font','arial','text',11,'labels',13,'box',1,'grid',1)
end


