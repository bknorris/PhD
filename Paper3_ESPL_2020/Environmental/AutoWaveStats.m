% This program automatically loads, pre-processes Aquadopp data, then runs
% a wavestats program to calculate wave statistics.
%
% Updates:
% 24/11/17: Instead of depth-averaging velocities immediately, I need the
% temporally averaged profiles to calculate u*c, the current shear
% velocity. I have also added two new fields to the runfile: 'ornt' for
% (up/down) and 'rd50' for the d50 of sediment used to calculate shear
% stresses. 
%
% BKN, UoW 2017
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear
%define working paths
dpath = '\Data\Aquadopp\';
ypath = {'Mekong_F2014';'Mekong_W2015'};
%load run file to tell program which files & VPs to load
rdir = 'd:\Projects\Mekong_W2015\DataAnalysis\Paper3\';
fid = fopen([rdir 'AutoAQDPrunfile.csv']);
rfile = textscan(fid,'%s%s%s%s%s%s%n','delimiter',',');
rexpt = rfile{1};rdate = rfile{2}; %experiment name, experiment date
rstart = rfile{3};rstop = rfile{4}; %start time, stop time (based on VPs)
rinst = rfile{5};rornt = rfile{6}; %instrument name in folder, instrument orientation
rd50 = rfile{7};                   %sediment d50 from Aaron's data
start = datenum(strcat(rdate,{' '},rstart),'dd-mm-yy HH:MM:SS');
stop = datenum(strcat(rdate,{' '},rstop),'dd-mm-yy HH:MM:SS');
stop = stop+datenum(0,0,0,0,3,0); %account for windowed indexing
for o = 1:2
    npath = ['g:\' ypath{o} dpath];
    folders = dir(npath);
    folders = {folders(3:end).name};
    for q = 1:length(folders)
        if any(strcmp(folders{q},rexpt))
            inid = strcmp(folders{q},rexpt);
            inid = find(inid);
            disp(['Loading files from ' npath folders{q} '\'])
        else
            continue
        end
        for j = 1:length(inid)
            disp(['Processing files from ' rdate{inid(j)}])
            afile = dir([npath folders{q} '\*' rinst{inid(j)} '*']);
            afile = {afile.name};
            disp(['Processing file: ' afile{1}])
            load([npath folders{q} '\' afile{:}])
            
            %%%Begin WaveStats routine%%%            
            disp(['Cropping data to times between ' datestr(start(inid(j)),'HH:MM:SS')...
                ' and ' datestr(stop(inid(j)),'HH:MM:SS')])
            lf = 1.2; %Cutoffs
            hf = 0.01;
            %Averaging interval
            win = 180; %seconds
            step = 1; %seconds
            tic
            zp = aqdp.metadata.HAB/1000;
            zuv = zp;
            time1 = aqdp.datenum;
            vid = find(time1 >= start(inid(j)) & time1 <= stop(inid(j)));
            time1 = time1(vid);
            U = nanmean(aqdp.u(vid),2);V = nanmean(aqdp.v(vid),2);W = nanmean(aqdp.w(vid),2);P = nanmean(aqdp.pressure(vid),2);
            x = (sin(aqdp.metadata.lat/57.29578))^2;
            fs = str2double(regexprep(aqdp.metadata.samprate,' Hz',''));
            avt = fs*step;
            nwin = fs*win;
            swin = fs*10; %10 second averaging window (to smooth)
            %DOF = round((nwin/swin)*2);
            
            U = cmgbridge(U,100,100,1000);V = cmgbridge(V,100,100,1000);
            W = cmgbridge(W,100,100,1000);P = cmgbridge(P,100,100,1000);
            U(isnan(U)) = 0;V(isnan(V)) = 0;W(isnan(W)) = 0;P(isnan(P)) = 0;
            
            %Save variables to structure
            wave.p = P;
            wave.u = U;
            wave.v = V;
            wave.w = W;
            wave.time = time1;
            wave.x = x;
            
            %Data Analysis
            disp('Calculating wave statistics...')
            nsamp = length(vid);
            ind = [1 avt:avt:nsamp];
            for ii = 1:length(ind)
                if abs(nsamp-ind(ii)) < nwin  %skip the last few indexes approaching the end of the t-s
                    continue
                else
                    idx = ind(ii):ind(ii)+nwin-1;
                end
                %Spectra of ADVs
                u = U(idx);
                v = V(idx);
                w = W(idx);
                p = P(idx)+zp;
                time = time1(ind(ii));
                
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
                
                %Wave parameters
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
                %Compute amplitude spectrum for direction/spreading
                Cpu = Cpu.*conj(Cpu);Cpv = Cpv.*conj(Cpv);
                Cuu = Cuu.*conj(Cuu);Cvv = Cvv.*conj(Cvv);
                Cuv = Cuv.*conj(Cuv);
                Dir=57.296*atan2(Cpu(lfc:hfc),Cpv(lfc:hfc));
                Dir=mod(Dir+180,360)-180;
                R2=((Cuu(lfc:hfc)-Cvv(lfc:hfc)).^2+4*Cuv(lfc:hfc).^2).^.5./(Cuu(lfc:hfc)+Cvv(lfc:hfc));
                Spread = 57.296*((1-R2)/2).^.5; %wave spreading
                omegar = sum(omega(lfc:hfc).*Guv(lfc:hfc)*df)./...
                    sum(Guv(lfc:hfc)*df); 
                ubr = sqrt(2*sum(Guv(lfc:hfc)*df)); %orbital wave velocity
                Tr = 2*pi/omegar; %orbital velocity period
                %Compute ripple parameters [Coco et al. 2007]
                
                
                %Save variables to structure
                wave.time2(ii) = time;
                wave.h(ii) = h;
                wave.Hs(ii) = Hs;
                wave.Tr(ii) = Tr;
                wave.F(ii,:) = F(lfc:hfc);
                wave.Dir(ii,:) = Dir;
                wave.Spread(ii,:) = Spread;
                wave.omegar(ii) = omegar;
                wave.ubr(ii) = ubr;
                wave.Uc(ii) = Uc;
                wave.Uwrms(ii) = Uwrms;
                wave.Wc(ii) = Wc;
                wave.Wwrms(ii) = Wwrms;
                wave.taub(ii) = taub;
            end
            disp(['Analysis completed in: ' num2str(toc/60) ' minutes'])
            iname = regexp(rinst{inid(j)},'[^_]+','match');
            sfpath = ['g:\' ypath{o} '\DataAnalysis\Paper3\WaveStats\' rdate{inid(j)} '\'];
            sfname = [rexpt{inid(j)} '_' rdate{inid(j)}(1:2) '_' iname{1}];
            disp(['Saving ' sfpath sfname '.mat'])
            save([sfpath sfname],'wave','-v7.3')
            clear wave aqdp
        end
    end
end
