%This program automatically loads, pre-processes Vectrino data, then
%calculates the roughness length scale of the bed (z0) and the bed shear
%stress from currents (tau_c). Script also loads the right wavestats file
%and combines the shear stress estimates using the Soulsby (1997) method.
clear
%define working paths
dpath = '\DataAnalysis\Paper3\';
ypath = {'Mekong_F2014';'Mekong_W2015'};
%load run file to tell program which files & VPs to load
rdir = 'd:\GradSchool\DataAnalysis\Paper3\ExperimentalDesign\';
fid = fopen([rdir 'AutoTaucRunfile_v2.csv']);
rfile = textscan(fid,'%s%s%s%n%n%s%n%n%n%n%n','delimiter',',');
rdate = rfile{1};rstart = rfile{2};%dates corresponding to folder names; start time for crop
rstop = rfile{3};rvp = rfile{4};  %stop time for crop; vp #
rhead = rfile{5};raqdp = rfile{6}; %VP heading for cross & along-shore rot;aquadopp file name to load
rsmpr= rfile{7};rhab = rfile{8}; %aqdp sample rate (Hz),vp height above bed (m)
rd50 = rfile{9};rvegd = rfile{10}; %d50, >0 if veg, 0 if no veg
rfd50 = rfile{11}; %force use of d50 for calculation of z0
start = datenum(strcat(rdate,{' '},rstart),'dd-mm-yy HH:MM:SS');
stop = datenum(strcat(rdate,{' '},rstop),'dd-mm-yy HH:MM:SS');
stop = stop+datenum(0,0,0,0,3,0); %account for windowed indexing
dfn = {'vpro1';'vpro2';'vpro3'};
tic
for i = 1:2
    %% Load the Data
    npath = ['e:\' ypath{i} dpath];
    folders = dir([npath '\VPs\']);
    folders = {folders(3:end).name};
    for ii = 1:length(folders)
        if any(strcmp(folders{ii},rdate))
            vpid = strcmp(folders{ii},rdate);
            vpid = find(vpid);
            disp(['Loading files from ' npath 'VPs\' folders{ii} '\'])
            disp(['Loading files from ' npath 'WaveStats\' folders{ii} '\'])
            disp(['Loading files from ' npath 'BottomTrack\' folders{ii} '\'])
            vfile = dir([npath 'VPs\' folders{ii} '\','*_Vels.mat']);
            bfile = dir([npath 'BottomTrack\' folders{ii} '\','*_bdtrace.mat']);
            bd = load([npath 'BottomTrack\' folders{ii} '\' bfile.name]);
        else
            continue
        end
        %Just load the variables you need, saves memory
        %% Process data for each Vectrino
        for j = 1:length(vpid)
            disp(['Processing ' dfn{rvp(vpid(j))}])
            disp(['Cropping data to times between ' datestr(start(vpid(j)),'HH:MM:SS')...
                ' and ' datestr(stop(vpid(j)),'HH:MM:SS')])
            data = load([npath 'VPs\' folders{ii} '\' vfile.name],dfn{rvp(vpid(j))});
            wfile = dir([npath 'WaveStats\' folders{ii} '\','*' raqdp{vpid(j)} '_v2.mat']);
            disp(['Loading ' npath 'WaveStats\' folders{ii} '\' wfile.name])
            load([npath 'WaveStats\' folders{ii} '\' wfile.name])
            %% Preliminary data processing
            vpt = data.(dfn{rvp(vpid(j))}).time;
            vid = find(vpt>=start(vpid(j))&vpt<=stop(vpid(j)));
            vpt = vpt(vid);
            vpx = data.(dfn{rvp(vpid(j))}).x(vid,:);
            vpy = data.(dfn{rvp(vpid(j))}).y(vid,:);
            bdt = bd.(dfn{rvp(vpid(j))}).time;
            bdist = bd.(dfn{rvp(vpid(j))}).bdist;
            rb = data.(dfn{rvp(vpid(j))}).rb;
            rho = 1010; %kg/m^3, estimated water density
            wvt = wave.time;
            [~,uid] = unique(bdt);
            bs = interp1(bdt(uid),bdist(uid),vpt);
            %Do the rotations from N-E to cross-along shore
            disp('Rotating to cross and along shore...')
            heading = rhead(vpid(j));th = heading*pi/180;
            R = [cos(th) -sin(th); sin(th) cos(th)];
            xr = zeros(size(vpx));yr = zeros(size(vpy));
            for k = 1:length(vpx)
                xx = vpx(k,:);
                yy = vpy(k,:);
                for kk = 1:length(xx)
                    rxy = [xx(kk) yy(kk)]*R;
                    xr(k,kk) = rxy(1);
                    yr(k,kk) = rxy(2);
                end
            end
            %% Calculate shear stresses
            %Average to the same length of time as in AutoWaveStats
            win = 180; %seconds
            step = 1; %seconds
            fs = 50;
            avt = fs*step;
            nwin = fs*win;
            disp('Calculating z0 and Tau_c...')
            nsamp = length(vid);
            ind = [1 avt:avt:nsamp];
            %Define parameters
            time = NaN(length(ind),1);
            z0 = NaN(length(ind),1);
            z = NaN(length(ind),35);
            phi = NaN(length(ind),1);
            umag = NaN(length(ind),35);
            xfit = NaN(length(ind),35);
            yfit = NaN(length(ind),35);
            discard = zeros(length(ind),1);
            for jj = 1:length(ind)
                if abs(nsamp-ind(jj)) < nwin  %skip the last few indexes approaching the end of the t-s
                    continue
                else
                    idx = ind(jj):ind(jj)+nwin-1;
                end
                %Time-average sections of time into profiles of velocity
                u = nanmean(sqrt(xr(idx,:).^2));
                vph = nanmean(bs(idx));
                zuv = rhab(vpid(j))-rb;
%                 zuv = vph-rb;
                zid = find(zuv<=0.005); %buffer zone of bed
                zuv(zid) = []; %remove below-bed values
                u(zid) = [];
                if rhab(vpid(j)) > 0.2 || rfd50(vpid(j)) %e.g., instruments too far above the bed
                    znot = (2.5*rd50(vpid(j))/30);
                    idd = NaN; maxrsq = 1; %signifies fits were not run
                    xs = NaN;
                    ys = NaN;
                    if jj == 1
                        flag = 'd50 used for z0';
                        disp(flag)
                    end
                else
                    if jj ==  1
                        flag = 'd50 estimated with log-fit';
                        disp(flag)
                    end
                    %Calculate z0 by fitting a line to the profile
                    if length(zuv) <= 2 || length(u) <= 2 %poor poly fits for 2 or less points!
                        znot = NaN;
                        idd = NaN; maxrsq = 1; %signifies fits were not run
                        xs = NaN;
                        ys = NaN;
                    else
                        %update: fit bottom n points while maintaining r^2 >
                        %0.8
                        xs = log(zuv);
                        ys = u;
                        minsamp = 5;
                        idd = linspace(length(zuv)-minsamp,1,length(zuv)-minsamp);
                        if isempty(idd)
                            discard(jj) = 1;
                            znot = NaN;
                            idd = NaN; maxrsq = 1;
                        else
                            rsqd = zeros(length(idd),1);
                            for q = 1:length(idd)
                                [~,~,~,~,stat] = regress(ys(idd(q):length(zuv))',...
                                    [ones(length(idd(q):length(zuv)),1) xs(idd(q):length(zuv))']);
                                rsqd(q) = stat(1);
                            end
                            if max(rsqd) < 0.8
                                znot = NaN;
                                idd = NaN; maxrsq = 1;
                            else
                                [~,maxrsq] = max(rsqd);
                                pf = polyfit(xs(idd(maxrsq):length(zuv)),ys(idd(maxrsq):length(zuv)),1);
                                znot = exp(-pf(2)/pf(1));
                            end
                            if znot > 1E-3 || znot < 1E-6 %discard z0 outside expected limits
                                znot = 2.5*rd50(vpid(j)); %nikuradse roughness length
                            end
                        end
                    end
                end
                %Calculate angle between currents and waves
                ux = nanmean(mean(xr(idx,:),2));
                vy = nanmean(mean(yr(idx,:),2));
                if ux == 0 || vy == 0
                    uwstat.phir = 0;
                else
                    uwstat = uvwstats(ux,vy);
                end
                %Save variables out to structure
                time(jj) = vpt(ind(jj));
                umag(jj,1:length(u)) = u;
                phi(jj) = uwstat.phir;
                z0(jj) = znot;
                z(jj,:) = vph-rb;
                xfit(jj,1:length(ys)) = ys;
                yfit(jj,1:length(xs)) = xs;
                clear uwstat
            end
            time(isnan(time))=[]; %remove trailing NaNs
            umag = umag(1:length(time),:);
            phi = phi(1:length(time));
            z0 = z0(1:length(time),:);
            if isnan(nanmean(z0)) 
                z0 = (2.5*rd50(vpid(j))/30);
                flag = 'd50 used for z0';
            else
                z0 = nanmean(z0);
                fprintf('z0 estimated as: %0.2d, d50 is: %0.2d\n',z0,rd50(vpid(j)))
            end  
            fprintf('There were %0.0f records discarded due to filtering protocol\n',nnz(discard))
            %Interpolate to the same length as wave.time2
            umag = interp1(time,umag,wave.time2);
            phi = interp1(time,phi,wave.time2);
            z = nanmean(z);
            z(z<0.005) = NaN;
            %Compute shear stresses
            %Update 31/10/2018: Use Yang & Nepf (2018) formulation for
            %veg shear stress in places where there was veg
            if rvegd(vpid(j)) > 0 %if there is veg
                disp('VP is in veg; using Y&N (2018) formula')
                Uc = nanmedian(umag,2);
                d = rvegd(vpid(j));
                nu = 1.05E-6; %kinematic viscosity
                Red = (Uc*d)/nu; %stem reynolds no.
                Cf = 1./((5.75*log10(2*wave.h/rd50(vpid(j)))).^2); %bed drag coefficient
                tauc = zeros(1,length(Uc));
                model1 = zeros(length(Uc),1); %for checking which model is used at the end of the run
                model2 = zeros(length(Uc),1);
                yang = 1;
                for jj = 1:length(Uc)
                    if Red(jj) < 4/Cf(jj)
                        tauc(jj) = (4*rho*nu*Uc(jj))/d;
                        model1(jj) = 1;
                    elseif Red(jj) >= 4/Cf(jj)
                        tauc(jj) = rho*Cf(jj)*Uc(jj)^2;
                        model2(jj) = 1;
                    end
                end
                %replace zeros with NaNs
                tauc(tauc == 0) = NaN;
            elseif rvegd(vpid(j)) == 0
                tauc = (rho.*((nanmedian(umag,2)*0.41)./(log(nanmedian(z)/z0))).^2)';
                yang = 0;
            end
            kw = 30*z0;
            azo = wave.ubr./(kw*wave.omegar);
            fw = exp((5.5.*azo.^-0.2)-6.3);
            tauw = 0.5*rho.*fw.*(wave.ubr.^2);
            taum = tauc.*(1+1.2*(tauw./(tauc+tauw)).^(3/2));
            tmax = sqrt((taum+tauw.*cos(phi)).^2+(tauw.*sin(phi)).^2);

            %% Save variables out to structure
            dat.time = wave.time2';
            dat.umag = umag;
            dat.z0 = z0;
            dat.z = z;
            dat.phi = phi';
            dat.tauc = tauc';
            dat.tauw = tauw';
            dat.tmax = tmax';
            dat.xfit = xfit;
            dat.yfit = yfit;
            dat.msg = flag;
            if yang == 1
                dat.model1 = sprintf('Percent of samples using (4*rho*nu*Uc)/d: %0.1f%%',...
                    (nnz(model1)/(nnz(model1)+nnz(model2)))*100); %percent of model 1
                dat.model2 = sprintf('Percent of samples using rho*Cf*Uc^2: %0.1f%%',...
                    (nnz(model2)/(nnz(model1)+nnz(model2)))*100); %percent of model 2
            end
            %Save file to folder
            iname = dfn{rvp(vpid(j))};
            expt = regexp(wfile.name,'[^_]+','match');
            sfpath = ['e:\' ypath{i} '\DataAnalysis\Paper3\BedStress\' rdate{vpid(j)} '\'];
            sfname = [expt{1} '_' rdate{vpid(j)}(1:2) '_' iname '_v2'];
            fprintf('Saving %s%s.mat\n\n',sfpath,sfname)
            save([sfpath sfname],'dat','-v7.3')
            clear wave dat data
        end
    end
end
disp(['Analysis completed in: ' num2str(toc/60) ' minutes'])
