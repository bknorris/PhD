% This program automatically loads, crops and runs the structure function
% on vertical vectrino velocities to estimate the dissipation rate of
% turbulence.
%
%
% BKN, UoW 2017
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear
%define working paths
dpath = '\DataAnalysis\Paper3\VPs\';
ypath = {'Mekong_F2014';'Mekong_W2015'};
%load run file to tell program which files & VPs to load
rdir = 'e:\MemoryStick\GradSchool\DataAnalysis\Paper3\ExperimentalDesign\';
fid = fopen([rdir 'AutoTKErunfile_v3.csv']);
rfile = textscan(fid,'%s%s%s%s%n%n','delimiter',',');
rexpt = rfile{1};rdate = rfile{2}; %experiment name, experiment date
rstart = rfile{3};rstop = rfile{4}; %start time, stop time (based on VPs)
rinst = rfile{5}; %instrument name in folder, instrument orientation
start = datenum(strcat(rdate,{' '},rstart),'dd-mm-yy HH:MM:SS');
stop = datenum(strcat(rdate,{' '},rstop),'dd-mm-yy HH:MM:SS');
stop = stop+datenum(0,0,0,0,3,0); %account for windowed indexing
disp(['Start time of run: ' datestr(now,'HH:MM:SS')])
for i = 1:2
    npath = ['e:\' ypath{i} dpath];
    folders = dir(npath);
    folders = {folders(3:end).name};
    for ii = 1:length(folders)
        if any(strcmp(folders{ii},rdate))
            inid = strcmp(folders{ii},rdate);
            inid = find(inid);
            disp(['Loading files from ' npath folders{ii} '\'])
        else
            continue
        end
        vpfile = dir([npath folders{ii} '\','*_Vels.mat']);
        file = matfile([npath folders{ii} '\' vpfile.name]);
        fn = whos('-file',[npath folders{ii} '\' vpfile.name]);
        fn = {fn.name};
        tic
        %% Run structure function code
        for j = 1:length(inid)
            disp(['Processing: ' fn{rinst(inid(j))}])
            data = file.(fn{rinst(inid(j))});
            %Basic structure function windowing settings
            fs = 50;
            window = 180; %180 second (3 min) averaging interval
            step = 1; %1 second step
            avt = step*fs; %samples/step
            nwin = window*fs; %samples/window
            cellsize = 0.9693;
            r = (cellsize/1000)/cosd(30); %converts vertical beam distance to along beam distance.                                                                    %vertical velocities are not converted.
            lags = 5;
            %Initialize structures for dataa storage
            disp('Running windowed structure function on vertical velocity bins')
            TKE = struct();
            nsamp = length(data.time);
            ind = [1 avt:avt:nsamp];
            [~,m] = size(data.x);
            maxbin = m-lags;       %TKE estimates are maximally length(maxbin) due to lag #
            for jj = 1:length(ind) %loop through time, windowed to the settings
                if abs(nsamp-ind(jj)) < nwin %skip the last few indexes approaching the end of the t-s
                    continue
                else
                    idx = ind(jj):ind(jj)+nwin-1;
                    DATA.z1 = data.z1(idx,:); %we are now using vertical beams only so Cv2 = 2.1.
                    DATA.z2 = data.z2(idx,:);
                    dfn = fieldnames(DATA);
                end
                if jj == 1
                    count = 1;
                elseif jj > 1
                    count = count+lags;
                end
                %Extract time of average
                time = data.time(ind(jj));
                for k = 1:length(dfn)                                   %loop through beams
                    %% Define variables:
                    itt = 1;                                            %iteration number
                    E = zeros(maxbin,1);
                    D = zeros(maxbin,lags);
                    pf = zeros(maxbin,2);
                    pv = zeros(maxbin,lags);
                    pf2 = zeros(maxbin,2);
                    pv2 = zeros(maxbin,lags);
                    pval = zeros(maxbin,1);
                    rsq = zeros(maxbin,1);
                    rsq2 = zeros(maxbin,1);
                    bfit = zeros(maxbin,1);
                    R = zeros(maxbin,lags);
                    rr = zeros(1,lags);
                    N = zeros(1,maxbin);
                    epsilon = zeros(1,maxbin);
                    epserr = zeros(1,maxbin);
                    beam = detrend(DATA.(dfn{k}),'constant');
                    %% Loop in Depth
                    while itt <= maxbin
                        d = zeros(length(idx),lags);
                        idl = itt:itt+lags;
                        idt = 1:length(idx);
                        for kk = 1:lags %compute velocity differences (vel(z,1) - vel(z,2:lags))^2
                            d(:,kk) = (beam(idt,itt)-(beam(idt,idl(kk+1)))).^2;
                            D(itt,kk) = nanmean(d(:,kk));
                            rr(:,kk) = r*kk;
                        end
                        R(itt,:) = rr.^(2/3); %along beam distance r^2/3
                        
                        [pf(itt,:),S] = polyfit(R(itt,:),D(itt,:),1);  %linear regression of D along r^2/3
                        A = pf(itt,1);
                        N(itt) = pf(itt,2);
                        epsilon(itt) = (A/2)^(3/2);  %units of W/m^3
                        itt = itt+1;
                    end
                    E(1:length(epsilon),1) = epsilon'; %TKE dissipation rate
                    E(abs(imag(E))>0) = NaN; %filter non-real numbers from E estimates, Noise
                    E(abs(imag(N'))>0) = NaN;
                    %% Save variables:
                    if k == 1
                        dat.(fn{rinst(inid(j))}).time(:,jj) = time;
                    end
                    dat.(fn{rinst(inid(j))}).(dfn{k}).E(1:maxbin,jj) = E;
                    dat.(fn{rinst(inid(j))}).(dfn{k}).N(1:maxbin,jj) = N;
                    dat.(fn{rinst(inid(j))}).(dfn{k}).lags = lags;
                    dat.(fn{rinst(inid(j))}).(dfn{k}).intv = [num2str(step) ' second step'];
                end
            end
            
        end
        disp(['Analysis completed in: ' num2str(toc/60) ' minutes'])
        iname = fn{rinst(inid(j))};
        sfpath = ['e:\' ypath{i} '\DataAnalysis\Paper3\Turbulence\' rdate{inid(j)} '\'];
        sfname = [rexpt{inid(j)} '_' rdate{inid(j)}(1:2) '_TKE_v3'];
        disp(['Saving ' sfpath sfname '.mat'])
        save([sfpath sfname],'-struct','dat','-v7.3')
    end
end