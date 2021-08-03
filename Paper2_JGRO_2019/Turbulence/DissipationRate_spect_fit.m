%Calculate dissipation rate from spectrum. This method employs the
%Trowbridge & Elgar (2001) method. Code used for thesis reviews.
%This is version 1 of this script, based on ReynoldsStressCFMethod_V4.mat
%
%
% Written by Benjamin K. Norris, University of Waikato 2019
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear
datdir = 'e:\Mekong_W2015\DataAnalysis\Paper3\VPs\';
files = {'7March2015_Vels.mat';'8March2015_Vels.mat';'10March2015_Vels.mat';...
    '14March2015a_Vels.mat'};
TKE = struct();
days = {'day1';'day2';'day3';'day4'};
%These are from PCA of the velocities to define the mean current direction
heading = [55 240 55 240]; %headings are in order! (HTA1 - 3, VTA)
lcutoff = [1.4 1.35 1.04 1.09]; %from wave cutoff freq (see CutoffFreqs)
ucutoff = [8.3 8.3 8.3 5];
for i = 2
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
    for ii = 1:3
        disp(['Running analysis on ' fn{ii}])
        dat = mfile.(fn{ii});
        lcf = lcutoff(i);
        ucf = ucutoff(i);
        disp(['Applying wave cutoff: ' sprintf('%0.2f',lcf) 'Hz'])
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
        disp(['Loading & rotations completed at: ' datestr(now,'HH:MM:SS')])
        x = b1;y = b2;z1 = b3;z2 = b4;
        
        %loop in time
        fs = 50;
        avt = fs*120; %2 minute step
        nwin = fs*600; %10 minute window (25% overlap)
        swin = fs*40; %40 second averaging window (to smooth spectra)
        nsamp = length(time);
        ind = [1 avt:avt:nsamp];
        for j = 1:length(ind)
            if abs(nsamp-ind(j)) < nwin  %skip the last few indexes approaching the end of the t-s
                continue
            else
                idx = ind(j):ind(j)+nwin-1;
            end
            time2 = time(idx(1));
            %initialize variables
            Euv = zeros(1,35);
            Eww = zeros(1,35);
            Euv_sd = zeros(1,35);
            Eww_sd = zeros(1,35);
            Zeta = zeros(1,35);
            for jj = 1:35
                X = x(idx,jj);
                Y = y(idx,jj);
                Z1 = z1(idx,jj);
                Z2 = z2(idx,jj);
                Z = (Z1+Z2)./2;
                X = detrend(X);Y = detrend(Y);Z1 = detrend(Z1);Z2 = detrend(Z2);
                
                %spectral settings
                nfft = nwin;
                noverlap = nfft/2;
                window = hamming(nfft);
                nf = (nfft/2)+1;
                [U,f]=pwelch(X,window,noverlap,nfft,fs); U = U(2:nf,1);
                [V,~]=pwelch(Y,window,noverlap,nfft,fs); V = V(2:nf,1);
                [W,~]=pwelch(Z,window,noverlap,nfft,fs); W = W(2:nf,1);
                f = f(2:nf,1);
                %use ubr instead of mean speed in this case
                Guv = U + V;
                range = find(f>0.05&f<lcf); %only in waveband
                df = f(3)-f(2);
                S = sqrt(2*sum(Guv(range)*df)); %wave orbital velocity in waveband
                
                %dissipation rate estimates
                twopi = 2*pi;

                % Frequency range that is used for inertial-range fit
                fm = find(f>lcf&f<ucf);
                fhn = find(f>ucf); % assume noise at f > 8 Hz
                wnoise  = min(W(fhn,1));
                N = wnoise*ones(size(W));
                W = W-N; %subtract off noise
                p = polyfit(log(f(fm)), log(W(fm)),1);
                zeta = p(1); %slope of intertial subrange fit
                Aww = p(2);
                %eq for epsilon from Tennekes & Lumeley, 1972
%                 B = twopi*((55/36)^(3/2)); %NB: this value is nearly 10
                B = 8; %e.g., Hay 2008
                eww = mean(B*((W(fm).^(3/2))/S).*(f(fm).^(5/2)));
                eww_sd = std(B*((W(fm).^(3/2))/S).*(f(fm).^(5/2)));
                
                %save data
                Eww(jj) = eww;
                Eww_sd(jj) = eww_sd;
                Zeta(jj) = zeta;
            end
            TKE.(date).(fn{ii}).time(j) = time2;
            TKE.(date).(fn{ii}).eww(:,j) = Eww;
            TKE.(date).(fn{ii}).eww_sd(:,j) = Eww_sd;
            TKE.(date).(fn{ii}).eww_slope(:,j) = Zeta;
        end
        clear dat
    end
    disp(['Processing time ' num2str(toc/60) ' minutes'])
end
disp(['Files completed at: ' datestr(now)])
fprintf('End of Run\n\n')
savedatdir = 'e:\Mekong_W2015\DataAnalysis\Paper2\Turbulence\';
fname = 'Dissip_spectral_lcf_ucf';
save([savedatdir fname],'TKE','-v7.3')
