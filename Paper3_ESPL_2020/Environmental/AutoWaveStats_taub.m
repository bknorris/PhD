% This program automatically loads, pre-processes Aquadopp data, then runs
% a wavestats program to calculate wave statistics.
%
% Updates:
% 24/11/17: Instead of depth-averaging velocities immediately, I need the
% temporally averaged profiles to calculate u*c, the current shear
% velocity. I have also added three new fields to the runfile: 'head' for
% the instrument heading, 'ornt' for (up/down) and 'rd50' for the d50 of 
% sediment used to calculate shear stresses. 
%
% BKN, UoW 2017
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear
%define working paths
dpath = '\Data\Aquadopp\';
ypath = {'Mekong_F2014';'Mekong_W2015'};
%load run file to tell program which files & VPs to load
rdir = 'd:\MemoryStick\GradSchool\DataAnalysis\Paper3\ExperimentalDesign\';
fid = fopen([rdir 'AutoAQDPrunfile.csv']);
rfile = textscan(fid,'%s%s%s%s%s%s%n%n%n','delimiter',',');
rexpt = rfile{1};rdate = rfile{2}; %experiment name, experiment date
rstart = rfile{3};rstop = rfile{4}; %start time, stop time (based on VPs)
rinst = rfile{5};rornt = rfile{6}; %instrument name in folder, instrument orientation
rhead = rfile{7};rd50 = rfile{8};  %transect heading (deg), sediment d50 from Aaron's data
rvegd = rfile{9};                  %0 if no veg; >0 if veg was present
start = datenum(strcat(rdate,{' '},rstart),'dd-mm-yy HH:MM:SS');
stop = datenum(strcat(rdate,{' '},rstop),'dd-mm-yy HH:MM:SS');
stop = stop+datenum(0,0,0,0,3,0); %account for windowed indexing
for i = 1:2
    npath = ['d:\' ypath{i} dpath];
    folders = dir(npath);
    folders = {folders(3:end).name};
    for ii = 1:length(folders)
        if any(strcmp(folders{ii},rexpt))
            inid = strcmp(folders{ii},rexpt);
            inid = find(inid);
            disp(['Loading files from ' npath folders{ii} '\'])
        else
            continue
        end
        for j = 1:length(inid)
            disp(['Processing files from ' rdate{inid(j)}])
            afile = dir([npath folders{ii} '\*' rinst{inid(j)} '*']);
            afile = {afile.name};
            disp(['Processing file: ' afile{1}])
            load([npath folders{ii} '\' afile{1}])
            
            %% Begin Wave Statistics Routine           
            disp(['Cropping data to times between ' datestr(start(inid(j)),'HH:MM:SS')...
                ' and ' datestr(stop(inid(j)),'HH:MM:SS')])
            tic
            %Extract key information from the Aquadopp metadata file
            x = (sin(aqdp.metadata.lat/57.29578))^2;
            fs = str2double(regexprep(aqdp.metadata.samprate,' Hz',''));
            zp = aqdp.metadata.HAB/1000; %height of pressure sensor
            blankdist = aqdp.metadata.blankdist;
            ncells = length(aqdp.rangebins);
            cellsize = aqdp.metadata.cellsize;
            %WaveStats parameters
            lf = 1.2; %Cutoffs
            hf = 0.01;
            win = 180; %seconds
            step = 1; %seconds
            avt = fs*step;
            nwin = fs*win;
            swin = fs*10; %10 second averaging window (to smooth)
            if j==1,DOF = round((nwin/swin)*2);
            disp(['Spectra averaged with: ' num2str(DOF) ' Degrees of freedom']),end
            %height of velocity measurements
            if strcmp(rornt{inid(j)},'up') 
                zuv = zp+blankdist+linspace(0,ncells*cellsize,ncells);
            elseif strcmp(rornt{inid(j)},'down')
                zuv = zp-blankdist-linspace(0,ncells*cellsize,ncells);
            end
            time1 = aqdp.datenum;
            vid = find(time1 >= start(inid(j)) & time1 <= stop(inid(j)));
            time1 = time1(vid);
            U = aqdp.u(vid,:);V = aqdp.v(vid,:);W = aqdp.w(vid,:);
            P = aqdp.pressure(vid,:);T = aqdp.temperature(vid,:);
            %bridge gaps in t-s with cmgbridge
            U = cmgbridge(U,100,100,1000);V = cmgbridge(V,100,100,1000);
            W = cmgbridge(W,100,100,1000);P = cmgbridge(P,100,100,1000);
            T = cmgbridge(T,100,100,1000);
            U(isnan(U)) = 0;V(isnan(V)) = 0;W(isnan(W)) = 0;
            P(isnan(P)) = 0;T(isnan(T)) = 0;
            %rotate to along & cross-shore
            heading = -360+rhead(inid(j));th = heading*pi/180; %CCW rotation!
            R = [cos(th) -sin(th); sin(th) cos(th)];
            Ur = zeros(size(U));Vr = zeros(size(V));
            for n = 1:ncells
                xy = [U(:,n) V(:,n)];
                for o = 1:length(xy);
                    rxy = [xy(o,1) xy(o,2)]*R;
                    Ur(o,n) = rxy(1);Vr(o,n) = rxy(2);
                end
            end
            %Save out original variables to structure
            wave.time = time1;
            wave.p = P;
            wave.u = Ur;
            wave.v = Vr;
            wave.w = W;
            
            %% Run the wave statistics analysis
            disp('Calculating wave statistics...')
            nsamp = length(vid);
            ind = [1 avt:avt:nsamp];
            for jj = 1:length(ind)
                if abs(nsamp-ind(jj)) < nwin  %skip the last few indexes approaching the end of the t-s
                    continue
                else
                    idx = ind(jj):ind(jj)+nwin-1;
                end
                %Spectra of ADVs
                u = nanmean(Ur(idx,:),2);
                v = nanmean(Vr(idx,:),2);
                w = nanmean(W(idx,:),2);
                p = P(idx)+zp;
                t = T(idx);
                time = time1(ind(jj));
                if strcmp(rornt{inid(j)},'up')
                    u0 = Ur(idx,1); %to estimate Cd (below)
                elseif strcmp(rornt{inid(j)},'down')
                    u0 = Ur(idx,end);
                end
                %Calculate RMS velocities (Luhar et al. 2013)
                Ec = (1/nwin)*sum(u);Nc = (1/nwin)*sum(v);
                Ewrms = sqrt((1/nwin)*sum((u-Ec).^2));
                Nwrms = sqrt((1/nwin)*sum((v-Nc).^2));
                Wc = (1/nwin)*sum(w);Wwrms = sqrt((1/nwin)*sum((w-Wc).^2));
                Uc = sqrt(Ec^2+Nc^2);Uwrms = sqrt(Ewrms^2+Nwrms^2);
                %     Uw = sqrt(2)*Uwrms;
                %     Ww = sqrt(2)*Wwrms;
                
                g = zeros(length(p),1);h = zeros(length(p),1);
                for l = 1:length(p)
                    g(l,:) = 9.780318*(1.0+(5.2788E-3+2.36E-5*x)*x)+1.092E-6*p(l,:);
                    h(l,:) = ((((-1.82E-15*p(l,:)+2.279E-10)*p(l,:)-2.2512E-5)*p(l,:)+9.72659)*p(l,:))/g(l,:);
                end
                rho = SeaDensity(15,mean(t),mean(p)+10);
                h = mean(h);
                g = mean(g);
                p = detrend(p);
                u = detrend(u);
                v = detrend(v);
                w = detrend(w);
                %Sig wave height should be calculated with non-directional
                %spectra, i.e. from the pressure t-s
                [Cpp,F] = pwelch(p,hanning(swin),swin*0.7,swin,fs);
                [Cuu,~] = pwelch(u,hanning(swin),swin*0.7,swin,fs);
                [Cvv,~] = pwelch(v,hanning(swin),swin*0.7,swin,fs);
                [Cww,~] = pwelch(w,hanning(swin),swin*0.7,swin,fs);
                [Cpu,~] = cpsd(p,u,hanning(swin),swin*0.7,swin,fs);
                [Cpv,~] = cpsd(p,v,hanning(swin),swin*0.7,swin,fs);
                [Cuv,~] = cpsd(u,v,hanning(swin),swin*0.7,swin,fs);
                Guv = Cuu+Cvv;
                
                %% Wave parameters
                lfc = find(F >= hf,1,'first');hfc = find(F <= lf,1,'last');
                df = F(3)-F(2);
                omega = 2*pi.*F;
                k = qkhf(omega,h)./h;
                kh = k*h;
                kz = k*median(zuv); %use median for depth-averaged velocities
                attn = cosh(kz)./cosh(kh);
                attn(attn<0.2) = 0.2;
                Spp = Cpp./(attn.^2);           %surface elevation spectrum
                m0 = sum(Spp(lfc:hfc)*df);
                Hs = 4*sqrt(m0);                %sig. wave height
                %Compute amplitude spectrum for direction/spreading
                Cpu = Cpu.*conj(Cpu);Cpv = Cpv.*conj(Cpv);
                Cuu = Cuu.*conj(Cuu);Cvv = Cvv.*conj(Cvv);
                Cuv = Cuv.*conj(Cuv);
                Dir=57.296*atan2(Cpu(lfc:hfc),Cpv(lfc:hfc));
                Dir=mod(Dir+180,360)-180;
                R2=((Cuu(lfc:hfc)-Cvv(lfc:hfc)).^2+4*Cuv(lfc:hfc).^2).^.5./(Cuu(lfc:hfc)+Cvv(lfc:hfc));
                Spread = 57.296*((1-R2)/2).^.5; %wave spreading
                omegar = sum(omega(lfc:hfc).*Guv(lfc:hfc)*df)./...
                    sum(Guv(lfc:hfc)*df); %orbital wave radian frequency
                ubr = sqrt(2*sum(Guv(lfc:hfc)*df)); %orbital wave velocity                
                Tr = 2*pi/omegar; %orbital velocity period
                %% Compute bed shear stress (Soulsby, 1997) 
                %Waves
                ks = 2.5*rd50(inid(j)); %Nikuradse roughness length (Soulsby, 1997)
                fw = exp(5.5*((ubr/ks*omegar)^-0.2)-6.3); %wave friction factor (Nielsen, 1992)
                tauw = 0.5*rho*fw*(ubr^2); %bed shear stress from waves
                %Currents 
                %Update 31/10/2018: Use Yang & Nepf (2018) formulation for
                %veg shear stress in places where there was veg
                if rvegd(inid(j)) > 0 %if there is veg
                    d = rvegd(inid(j));
                    nu = 1.05E-6; %kinematic viscosity
                    Red = (Uc*d)/nu; %stem reynolds no.
                    Cf = 1/((5.75*log10(2*h/rd50(inid(j))))^2); %bed drag coefficient
                    if Red < 4/Cf
                        tauc = (4*rho*nu*Uc)/d;
                        disp('no')
                    elseif Red >= 4/Cf
                        tauc = rho*Cf*Uc^2;
                        disp('yes')
                    end
                elseif rvegd(inid(j)) == 0
                    ustar = (Uc*0.4)/(log(median(zuv)/ks));
                    tauc = rho*(ustar^2);
                    Cf = (ustar^2)/mean(u0)^2; %Lacy et al. (2011)
                end
                %Combined [note, ignores influence of vegetative drag]
                uwstat = uvwstats(u,v); 
                phi = uwstat.phir;
                taum = tauc*(1+1.2*(tauw/(tauc+tauw))^3.2);
                tmax = ((taum+tauw*cos(phi))^2+(tauw*sin(phi))^2)^(1/2);
                %% Save variables to structure
                wave.time2(jj) = time; 
                wave.h(jj) = h; %depth [m]
                wave.Hs(jj) = Hs; %sig. wave height [m]
                wave.rho(jj) = rho; %water density [kg/m3]
                wave.Cd(jj) = Cf; %estimated drag coefficient [-]
                wave.Tr(jj) = Tr; %orbital velocity period [s]
                wave.F(jj,:) = F(lfc:hfc); %frequency band between cutoff freqs.
                wave.Dir(jj,:) = Dir; %wave direction [deg]
                wave.Spread(jj,:) = Spread; %wave spreading [deg]
                wave.omegar(jj) = omegar; %wave radian frequency [rad/s]
                wave.fw(jj) = fw; %wave friction factor [-]
                wave.ubr(jj) = ubr; %wave orbital velocity in freq. band [m/s]
                wave.Uc(jj) = Uc; %velocity magnitude of currents, horizontal [m/s]
                wave.Uwrms(jj) = Uwrms; %velocity magnitude of waves, horizontal [m/s]
                wave.Wc(jj) = Wc; %velocity magnitude of currents, vertical [m/s]
                wave.Wwrms(jj) = Wwrms; %velocity magnitude of waves, vertical [m/s]
                wave.tauc(jj) = tauc; %shear stress due to currents [kg/ms2]
                wave.tauw(jj) = tauw; %shear stress due to waves [kg/ms2]
                wave.tmax(jj) = tmax; %maximum combined current and wave shear stress [kg/ms2]
            end
            disp(['Analysis completed in: ' num2str(toc/60) ' minutes'])
            iname = regexp(rinst{inid(j)},'[^_]+','match');
            sfpath = ['d:\' ypath{i} '\DataAnalysis\Paper3\WaveStats\' rdate{inid(j)} '\'];
            sfname = [rexpt{inid(j)} '_' rdate{inid(j)}(1:2) '_' iname{1}];
            disp(['Saving ' sfpath sfname '.mat'])
            save([sfpath sfname],'wave','-v7.3')
            clear wave aqdp
        end
    end
end
