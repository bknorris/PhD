%Automate velocity averaging for VPs using the same windowing as the
%Structure Function
clear
maindir = 'd:\Projects\Mekong_W2015\DataAnalysis\Paper1\';
savedatdir = 'd:\Projects\Mekong_W2015\DataAnalysis\Paper2\AveragedVelocities\';
files = dir([maindir '*.mat']);files = {files.name};
toprocess = 10:16;
year = '2015'; %'2014' or '2015'
%to set the headings, load the transects & rotation data from a file
fid = fopen([maindir 'transects.txt']);
trans = textscan(fid,'%n','delimiter',',');
trans = trans{1};
fid = fopen([maindir 'VPheadings.txt']);
head = textscan(fid,'%n','delimiter',',');
head = head{1};
if strcmp(year,'2014')
    fid = fopen([maindir 'horizontal.txt']);
    horiz = textscan(fid,'%n','delimiter',',');
    horiz = horiz{1};horiz = reshape(horiz,[3,8]);
    trans = reshape(trans,[3,8]);head = reshape(head,[3,8]);
elseif strcmp(year,'2015')
    trans = reshape(trans,[3,16]);head = reshape(head,[3,16]);
    horiz = zeros(3,16);
end
%calculate headings
headings = (head+360)-trans;
%%%%
for k = 1:length(toprocess)
    disp(['Loading ' files{toprocess(k)}])
    load([maindir files{toprocess(k)}])
    
    fn = fieldnames(dat);
    %%%%
    window = 60; %60 second averaging interval
    step = 30; %30 second step
    fs = 50;
    avt = step*fs; %samples/step
    nwin = window*fs; %samples/window
    Vels = struct();                                                            %initialize structures for data storage
    Avgs = struct();
    %%%%
    [tp,~] = size(fn);
    if tp == 1
        g = 1;
    else
        g = 1:3;
    end
    for i = g                                                                 %loop through instruments VP1-VP3
        disp(['Analysing ' fn{i}])
        nsamp = length(dat.(fn{i}).time);
        ind = [1 avt:avt:nsamp];
        [~,m] = size(dat.(fn{i}).beam1);
        %%%%
        for ii = 1:length(ind)                                                  %loop through time, windowed to the settings
            if abs(nsamp-ind(ii)) < nwin                                        %skip the last few indexes approaching the end of the t-s
                continue
            else
                idx = ind(ii):ind(ii)+nwin-1;
                x = dat.(fn{i}).x(idx,:);
                y = dat.(fn{i}).y(idx,:);
                z1 = dat.(fn{i}).z1(idx,:);
                z2 = dat.(fn{i}).z2(idx,:);
            end
            if strcmp(year,'2014') && horiz(i,toprocess(k)) == 1               %IF the instruments are horizontal, swap z and x beams, invert y and x beams for rotations
                if ii == 1
                    disp('User reports the instrument is HORIZONTAL')
                    disp('Swapping z and x beams')
                end
                z = x;
                x = (z1+z2)./2;
                x = x*-1;y = y*-1;
            end
            heading = headings(i,toprocess(k));
            if ii == 1
                disp(['Rotating instrument to specified heading: ' num2str(heading) sprintf('%c', char(176))])
            end
            %rotate to cross-shore and along-shore
            rot = heading*pi/180;
            x = x.*(ones(size(x))*cos(rot)) + ...
                y.*(ones(size(y))*sin(rot));
            y = -y.*(ones(size(y))*sin(rot)) + ...
                x.*(ones(size(x))*cos(rot));
            
            %extract time of average
            Vels.time(ii) = dat.(fn{i}).time(ind(ii));
            
            %calculate RMS velocities (Luhar et al. 2013)
            Uc = zeros(35,1);
            Uw = zeros(35,1);
            Umag = zeros(35,1);
            for itt = 1:35
                Ec = (1/length(y(:,itt)))*sum(y(:,itt));Nc = (1/length(x(:,itt)))*sum(x(:,itt));
                Ewrms = sqrt((1/length(y(:,itt)))*sum((y(:,itt)-Ec).^2));
                Nwrms = sqrt((1/length(x(:,itt)))*sum((x(:,itt)-Nc).^2));
                Uc(itt) = sqrt(Ec^2+Nc^2);Uwrms = sqrt(Ewrms^2+Nwrms^2);
                Uw(itt) = sqrt(2)*Uwrms;
                
                %calculate horizontal velocity magnitude
                Umag(itt) = nanmean(sqrt(x(:,itt).^2+y(:,itt).^2));
            end
            
            %compute |u|^3
            ucubed = (y.^2+x.^2).^(3/2);    %'x' points 'north', 'y' points 'west'
            
            %compute u', v' and w'
            up = 
            %calculate velocity gradient
            win = 5*fs;
            nn = length(x);
            ind2 = [1 fs:fs:nn];
            dudz = NaN(length(ind2),35);
            zh = linspace(0.001,0.035,35);
            for j = 1:length(ind2)
                if abs(nn-ind2(j)) < win
                    continue
                else
                    idx2 = ind2(j):ind2(j)+win-1;
                    u = nanmean(x(idx2,:));
                    if sum(u) == 0
                        dudz(j,:) = NaN;
                    else
                        ind3 = [1 3:2:35];
                        ind4 = 2:2:35;
                        for jj = 1:length(ind3)-1
                            b = polyfit(u(ind3(jj):ind3(jj+1)),zh(ind3(jj):ind3(jj+1)),1);
                            dudz(j,ind4(jj)) = b(1);
                        end
                    end
                end
            end
            if nnz(isnan(dudz))/nnz(dudz) < 1
               dudz = fixgaps(nanmean(dudz));
            else
                dudz = dudz(1,:);
            end
            
            %average in time per depth bin
            mx = nanmean(x)';
            my= nanmean(y)';
            mz1= nanmean(z1)';
            mz2 = nanmean(z2)';
            uc = nanmean(ucubed)';
            dUdZ = smooth(dudz,6);
            
            Vels.x(:,ii) = mx;
            Vels.y(:,ii) = my;
            Vels.z1(:,ii) = mz1;
            Vels.z2(:,ii) = mz2;
            Vels.Uc(:,ii) = Uc;
            Vels.Uw(:,ii) = Uw;
            Vels.Umag(:,ii) = Umag;
            Vels.ucubed(:,ii) = uc;
            Vels.dudz(:,ii) = dUdZ;
        end
        Avgs.(fn{i}) = Vels;
    end
    file = regexprep(files{toprocess(k)},'Vels.mat','');
    disp(['Saving ' file ' Averaged velocity file'])
    fname = [file 'Vave'];
    save([savedatdir fname],'Avgs','-v7.3')
end