clear
close all
datdir = 'd:\Projects\Mekong_W2015\DataAnalysis\Paper3\VPs\';
files = {'7March2015_Vels.mat';'8March2015_Vels.mat';'10March2015_Vels.mat';...
    '14March2015a_Vels.mat'};
RMS = struct();
days = {'day1';'day2';'day3';'day4'};
%These are from PCA of the velocities to define the mean current direction
heading = [55 240 55 260]; %headings are in order! (HTA1 - 3, VTA) was 55
locutoff = [1.4 1.35 1.04 1.09]; %from wave cutoff freq (see CutoffFreqs)
hicutoff = [8.3 8.3 8.3 5]; %high f cutoff
for i = 1:4
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
        lcf = locutoff(i);
        hcf = hicutoff(i);
        dfn = fieldnames(dat);
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
        %loop in time
        fs = 50;
        avt = fs*150; %2.5 minute step
        nwin = fs*600; %10 minute window (25% overlap)
        swin = fs*40; %10 second averaging window (to smooth spectra)
        nsamp = length(time);
        ind = [1 avt:avt:nsamp];
        for j = 1:length(ind)
            if abs(nsamp-ind(j)) < nwin  %skip the last few indexes approaching the end of the t-s
                continue
            else
                idx = ind(j):ind(j)+nwin-1;
            end
            time2 = time(idx(nwin/2));
            if i == 1 || i == 4
                bins = 1:5;
            else
                bins = 13:19;
            end
            X = x(idx,bins);
            Y = y(idx,bins);
            Z1 = z1(idx,bins);
            Z2 = z2(idx,bins);
            u = nanmean(detrend(X),2);
            v = nanmean(detrend(Y),2);
            w1 = nanmean(detrend(Z1),2);
            w2 = nanmean(detrend(Z2),2);
            w = (w1+w2)./2;
            %estimate rms velocities (sqrt of integrated spectrum)
            [Guu,f] = pwelch(u,swin,swin/2,[],fs);
            [Gvv,~] = pwelch(v,swin,swin/2,[],fs);
            [Gww,~] = pwelch(w,swin,swin/2,[],fs);
            range = find(f>lcf&f<hcf);
            urms = sqrt(trapz(Guu(range)));
            vrms = sqrt(trapz(Gvv(range))); 
            wrms = sqrt(trapz(Gww(range)));
            %save to file
            RMS.(date).(fn{ii}).time(j,:) = time2;
            RMS.(date).(fn{ii}).urms(j,:) = urms;
            RMS.(date).(fn{ii}).vrms(j,:) = vrms;
            RMS.(date).(fn{ii}).wrms(j,:) = wrms;
        end
    clear dat
    end
    disp(['Processing time ' num2str(toc/60) ' minutes'])
end
disp(['Files completed at: ' datestr(now)])
fprintf('End of Run\n\n')
savedatdir = 'd:\Projects\Mekong_W2015\DataAnalysis\Paper2\AveragedVelocities\';
save([savedatdir 'RMS_inertsubrng_vels'],'RMS','-v7.3')

            
            
