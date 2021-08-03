%Calculate average velocities for the turbulence statistics calculations.
clear
close all
datdir = 'd:\Projects\Mekong_W2015\DataAnalysis\Paper3\VPs\';
savedatdir = 'd:\Projects\Mekong_W2015\DataAnalysis\Paper2\AveragedVelocities\Two\';
files = {'7March2015_Vels.mat';'8March2015_Vels.mat';'10March2015_Vels.mat';...
    '14March2015a_Vels.mat'};
days = {'day1';'day2';'day3';'day4'};
heading = [55 240 55 260]; %headings are in order! (HTA1 - 3, VTA) was 55
Avgs = struct();
for i = 4
    tic
    disp(['Run started at: ' datestr(now)])
    disp(['Loading ' files{i}])
    fn = whos('-file',[datdir files{i}]);
    fn = {fn.name};
    mfile = matfile([datdir files{i}]);
    if length(fn) == 3
        g = 1:3;
    else
        g = 1;
    end
    for ii = g
        disp(['Running analysis on ' fn{ii}])
        dat = mfile.(fn{ii});
        time = dat.time;
        dfn = fieldnames(dat);
        %         for iii = 8:15
        %             dat = rmfield(dat,dfn{iii}); %clean some memory
        %         end
        if i < 3
            ind = 1:length(time)-1;
            date = days{i};
        elseif i == 3
            start = datenum(2015,03,10,15,14,00);
            stop = datenum(2015,03,10,16,30,00);
            ind = find(time >= start & time <= stop);
            time = time(ind);
            date = days{i};
        else
            start = datenum(2015,03,14,07,10,09);
            stop = datenum(2015,03,14,09,20,09);
            ind = find(time >= start & time <= stop);
            time = time(ind);
            date = days{4};
        end
        b1 = dat.x(ind,:);
        b2 = dat.y(ind,:);
        b3 = dat.z1(ind,:);
        b4 = dat.z2(ind,:);
        
        %Update: 12/02/18 Steve recommends rotating horizontal velocities to the
        %principal component direction
%         out = cmgpca(navg(nanmean(b1,2),1000),navg(nanmean(b2,2),1000),[],0);
%         fprintf('Rotating velocities to %0.2f degrees\n',out.mdir)
        fprintf('Rotating velocities to %0.2f degrees\n',heading(i))
        
%         th = out.mdir*pi/180;
        th = heading(i)*pi/180;
        R = [cos(th) -sin(th); sin(th) cos(th)];
        [m,n] = size(b1);
        rb1 = zeros(m,n);
        rb2 = zeros(m,n);
        for k = 1:35
            for kk = 1:length(b1)
                rxy = [b1(kk,k) b2(kk,k)]*R;
                rb1(kk,k) = rxy(1);
                rb2(kk,k) = rxy(2);
            end
        end
        b1 = rb1;b2 = rb2;clear rb1 rb2 rxy
        x = b1;y = b2;z1 = b3;z2 = b4;
        disp(['Loading & rotations completed at: ' datestr(now,'HH:MM:SS')])
        %loop in time
        fs = 50;
        avt = fs*150; %2.5 minute step
        nwin = fs*600; %10 minute window (25% overlap)
        nsamp = length(time);
        ind = [1 avt:avt:nsamp];
        Vels = struct();
        for j = 1:length(ind)
            if abs(nsamp-ind(j)) < nwin  %skip the last few indexes approaching the end of the t-s
                continue
            else
                idx = ind(j):ind(j)+nwin-1;
            end
            X = x(idx,:);
            Y = y(idx,:);
            Z1 = z1(idx,:);
            Z2 = z2(idx,:);
            %extract time of average
            Vels.time(j) = dat.time(ind(j));
            
            %calculate RMS velocities (Luhar et al. 2013)
            Uc = zeros(35,1);
            Uw = zeros(35,1);
            Umag = zeros(35,1);
            for itt = 1:35
                Ec = (1/length(Y(:,itt)))*sum(Y(:,itt));Nc = (1/length(X(:,itt)))*sum(X(:,itt));
                Ewrms = sqrt((1/length(Y(:,itt)))*sum((Y(:,itt)-Ec).^2));
                Nwrms = sqrt((1/length(X(:,itt)))*sum((X(:,itt)-Nc).^2));
                Uc(itt) = sqrt(Ec^2+Nc^2);Uwrms = sqrt(Ewrms^2+Nwrms^2);
                Uw(itt) = sqrt(2)*Uwrms;
                
                %calculate horizontal velocity magnitude
                Umag(itt) = nanmean(sqrt(X(:,itt).^2+Y(:,itt).^2));
            end
            
            %compute |u|^3
            ucubed = sqrt(Y.^2+X.^2).^(3/2);    %'x' points 'north', 'y' points 'west'

            %average in time per depth bin
            mx = nanmean(X)';
            my= nanmean(Y)';
            mz1= nanmean(Z1)';
            mz2 = nanmean(Z2)';
            uc = nanmean(ucubed)';

            
            Vels.x(:,j) = mx;
            Vels.y(:,j) = my;
            Vels.z1(:,j) = mz1;
            Vels.z2(:,j) = mz2;
            Vels.Uc(:,j) = Uc;
            Vels.Uw(:,j) = Uw;
            Vels.Umag(:,j) = Umag;
            Vels.ucubed(:,j) = uc;
            if j == 1
                Vels.intv = [num2str(avt/fs) ' second step'];
                Vels.win = [num2str(nwin/fs) ' second window'];
            end
        end
        Avgs.(fn{ii}) = Vels;
        disp([fn{ii} ' completed at: ' datestr(now,'HH:MM:SS')])
    end
    disp(['Processing time ' num2str(toc/60) ' minutes'])   
    file = regexprep(files{i},'Vels.mat','');
    disp(['Saving ' file ' Averaged velocity file'])
    fname = [file 'Vave'];
    save([savedatdir fname],'Avgs','-v7.3')
end
disp(['Files completed at: ' datestr(now)])
    
        
