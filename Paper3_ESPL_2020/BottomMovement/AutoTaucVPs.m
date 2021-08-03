%This program automatically loads, pre-processes Vectrino data, then
%calculates the roughness length scale of the bed (z0) and the bed shear
%stress from currents (tau_c). Script also loads the right wavestats file
%and combines the shear stress estimates using the Soulsby (1997) method.
clear
%define working paths
dpath = '\DataAnalysis\Paper3\';
ypath = {'Mekong_F2014';'Mekong_W2015'};
%load run file to tell program which files & VPs to load
rdir = 'e:\MemoryStick\GradSchool\DataAnalysis\Paper3\ExperimentalDesign\';
fid = fopen([rdir 'AutoTaucRunfile_v2.csv']);
rfile = textscan(fid,'%s%s%s%n%n%s%n%n%n','delimiter',',');
rdate = rfile{1};rstart = rfile{2};%dates corresponding to folder names; start time for crop
rstop = rfile{3};rvp = rfile{4};  %stop time for crop; vp #
rhead = rfile{5};raqdp = rfile{6}; %VP heading for cross & along-shore rot;aquadopp file name to load
rsmpr= rfile{7};rhab = rfile{8}; %aqdp sample rate (Hz),vp height above bed (m)
rd50 = rfile{9};                 %d50
start = datenum(strcat(rdate,{' '},rstart),'dd-mm-yy HH:MM:SS');
stop = datenum(strcat(rdate,{' '},rstop),'dd-mm-yy HH:MM:SS');
stop = stop+datenum(0,0,0,0,3,0); %account for windowed indexing
dfn = {'vpro1';'vpro2';'vpro3'};
tic
for i = 1:2
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
        for j = 1:length(vpid)
            disp(['Processing ' dfn{rvp(vpid(j))}])
            disp(['Cropping data to times between ' datestr(start(vpid(j)),'HH:MM:SS')...
                ' and ' datestr(stop(vpid(j)),'HH:MM:SS')])
            data = load([npath 'VPs\' folders{ii} '\' vfile.name],dfn{rvp(vpid(j))});
            wfile = dir([npath 'WaveStats\' folders{ii} '\','*' raqdp{vpid(j)} '.mat']);
            disp(['Loading ' npath 'WaveStats\' folders{ii} '\' wfile.name])
            load([npath 'WaveStats\' folders{ii} '\' wfile.name])
            %Crop to the start/stop times
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
            for jj = 1:length(ind)
                if abs(nsamp-ind(jj)) < nwin  %skip the last few indexes approaching the end of the t-s
                    continue
                else
                    idx = ind(jj):ind(jj)+nwin-1;
                end
                %Time-average sections of time into profiles of velocity
                u = nanmean(sqrt(xr(idx,:).^2));
                vph = nanmean(bs(idx));
                zuv = vph-rb;
                zid = find(zuv<=0.005); %buffer zone of bed
                zuv(zid) = []; %remove below-bed values
                u(zid) = [];
                %Calculate z0 by fitting a line to the profile
                if length(zuv) <= 2 || length(u) <= 2 %poor poly fits for 2 or less points!
                    znot = NaN;
                else
                    xs = log(zuv);
                    ys = u;
                    pf = polyfit(xs,ys,1);
                    znot = exp(-pf(2)/pf(1));
                end
                if znot > 0.1 %discard z0 greater than 10 cm
                    znot = NaN;
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
                clear uwstat
            end
            time(isnan(time))=[]; %remove trailing NaNs
            umag = umag(1:length(time),:);
            phi = phi(1:length(time));
            z0 = z0(1:length(time),:);
            if rhab(vpid(j)) > 0.01 
                z0 = rd50(vpid(j));
                disp('Using D50 for z0')
            end
            %Interpolate to the same length as wave.time2
            umag = interp1(time,umag,wave.time2);
            phi = interp1(time,phi,wave.time2);
            %Compute shear stresses
            z = nanmean(z);
            zuv = z;
            z(z<0.005) = NaN;
            tauc = (rho.*((nanmedian(umag,2)*0.4)./(log(nanmedian(z)/max(z0)))).^2)';
            ks = max(z0)*30;
            fw = exp((5.5.*(wave.ubr./(ks.*wave.omegar)).^-0.2)-6.3);
            tauw = 0.5*rho.*fw.*(wave.ubr.^2);
            taum = tauc.*(1+1.2*(tauw./(tauc+tauw)).^3.2);
            tmax = ((taum+tauw.*cos(phi)).^2+(tauw.*sin(phi)).^2).^(1/2);
            %Save variables out to structure
            dat.time = wave.time2;
            dat.umag = umag;
            dat.z0 = max(z0);
            dat.z = zuv;
            dat.phi = phi';
            dat.tauc = tauc;
            dat.tauw = tauw;
            dat.tmax = tmax;
            %Save file to folder
            iname = dfn{j};
            expt = regexp(wfile.name,'[^_]+','match');
            sfpath = ['d:\' ypath{i} '\DataAnalysis\Paper3\BedStress\' rdate{vpid(j)} '\'];
            sfname = [expt{1} '_' rdate{vpid(j)}(1:2) '_' iname];
            disp(['Saving ' sfpath sfname '.mat'])
            save([sfpath sfname],'dat','-v7.3')
            clear wave dat data
        end
    end
end
disp(['Analysis completed in: ' num2str(toc/60) ' minutes'])
