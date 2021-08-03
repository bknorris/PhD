%Calculate Reynolds stress by calculating the cospectra of velocities, then
%filtering by the turbulent frequency band previously defined by our other
%analyses. We are adapting the CF method of Gerbi et al., 2008 here, but
%are not focusing on model fits to derive Reynolds stress. Instead, we
%simply calculate cospectra, define the cutoff frequency, and integrate the
%cospectra above this frequency. This is version 2 of this script.
clear
close all
datdir = 'd:\Projects\Mekong_W2015\DataAnalysis\Paper1\';
files = {'HTA_1Vels.mat';'HTA_2Vels.mat';'HTA_4Vels.mat';...
    'VTA_2vp1Vels.mat';'VTA_2vp2Vels.mat';'VTA_2vp3Vels.mat'};
RS = struct();
days = {'day1';'day2';'day3';'day4'};
angle = 5; %rotation angle for testing limitations
orienttest = 0;
% cutoff = 4;
cutoff = [8.2 8.2 8.2;7 7 7;7 7 7;5.1 0 0 ;5.1 0 0;5.1 0 0]; %from wave cutoff freq (see WaveStatistics_v2)
tic
for i = 1:6
    disp(['Loading ' files{i}])
    load([datdir files{i}])
    fn = fieldnames(dat);
    if length(fn) == 3
        g = 1:3;
    else
        g = 1;
    end
    for ii = g
        disp(fn{ii})
        cf = cutoff(i,ii);
        disp(['Applying cutoff: ' sprintf('%0.1f',cf)])
        time = dat.(fn{ii}).time;
        dfn = fieldnames(dat.(fn{ii}));
        for k = 8:15
            dat.(fn{ii}) = rmfield(dat.(fn{ii}),dfn{k}); %clean some memory
        end
        if i < 3
            heading = 20;
            ind = 1:length(time)-1;
            date = days{i};
        elseif i == 3
            heading = 20;
            start = datenum(2015,03,10,15,14,00);
            stop = datenum(2015,03,10,16,30,00);
            ind = find(time >= start & time <= stop);
            time = time(ind);
            date = days{i};
        else
            heading = 96;
            start = datenum(2015,03,14,07,10,09);
            stop = datenum(2015,03,14,09,20,09);
            ind = find(time >= start & time <= stop);
            time = time(ind);
            date = days{4};
        end
        b1 = dat.(fn{ii}).x(ind,:);
        b2 = dat.(fn{ii}).y(ind,:);
        b3 = dat.(fn{ii}).z1(ind,:);
        b4 = dat.(fn{ii}).z2(ind,:);
        rot = (pi*heading)/180;
        b1 =  b1.*cos(rot)+b2.*sin(rot);
        b2 = -b1.*sin(rot)+b2.*cos(rot);
        
        if orienttest
            %perform an additional rotation to test if RS changes by
            %altering the angles by a small amount
            disp('User requests +/-5 deg pitch/roll adjustment to velocities')
            pp = 5*pi/180;
            rr = -5*pi/180;
            P = [cos(pp) -sin(pp)*sin(rr) -cos(rr)*sin(pp);...
                0             cos(rr)          -sin(rr);  ...
                sin(pp) sin(rr)*cos(pp)  cos(pp)*cos(rr)];
            [m,n] = size(b1);
            x = zeros(m,n);y = zeros(m,n);z1 = zeros(m,n);z2 = zeros(m,n);
            for j = 1:n
                V = [b1(:,j) b2(:,j) b3(:,j)];
                vels = V*P;
                x(:,j) = vels(:,1);
                y(:,j) = vels(:,2);
                z1(:,j) = vels(:,3);
                z2(:,j) = z2(:,j); %4-beam ADCP transforms ignore beam 4.
            end
        else
            x = b1;y = b2;z1 = b3;z2 = b4;
        end
        %loop in time
        fs = 50;
        avt = fs*150; %2.5 minute step
        nwin = fs*600; %10 minute window (25% overlap)
        swin = fs*10; %10 second averaging window (to smooth spectra)
        nsamp = length(time);
        ind = [1 avt:avt:nsamp];
        UW = zeros(1,35);
        VW = zeros(1,35);
        UWSTAR = zeros(1,35);
        VWSTAR = zeros(1,35);
        KU0 = zeros(1,35);
        KV0 = zeros(1,35);
        matlabpool(4)
        for j = 1:length(ind)
            if abs(nsamp-ind(j)) < nwin  %skip the last few indexes approaching the end of the t-s
                continue
            else
                idx = ind(j):ind(j)+nwin-1;
            end
            time2 = time(ind(j));
            parfor jj = 1:35
                X = x(idx,jj);
                Y = y(idx,jj);
                Z1 = z1(idx,jj);
                Z2 = z2(idx,jj);
                if min(X) == 0 && max(X) == 0 || min(Y) == 0 && max(Y) == 0 %skip bins with zero velocity
                    UW(jj) = 0;
                    VW(jj) = 0;
                    UWSTAR(jj) = 0;
                    VWSTAR(jj) = 0;
                    KU0(jj) = 0;
                    KV0(jj) = 0;
                else
                    u = detrend(X);
                    v = detrend(Y);
                    w1 = detrend(Z1);
                    w2 = detrend(Z2);
                    theta = 30;
                    
                    %calc along-beam velocities
                    u1 = -u*sind(theta)-w1*cosd(theta);
                    u2 = u*sind(theta)-w2*cosd(theta);
                    u3 = -v*sind(theta)-w1*cosd(theta);
                    u4 = v*sind(theta)-w2*cosd(theta);
                    
                    %compute velocity spectra
                    [Su1,fu] = pwelch(u1,swin,swin*0.7,nwin,fs,'onesided');
                    [Su2,~] = pwelch(u2,swin,swin*0.7,nwin,fs,'onesided');
                    [Su3,~] = pwelch(u3,swin,swin*0.7,nwin,fs,'onesided');
                    [Su4,~] = pwelch(u4,swin,swin*0.7,nwin,fs,'onesided');
                    
                    %calculate cospectra
                    COuw = (Su1-Su2)./(4*cosd(theta)*sind(theta));
                    COvw = (Su3-Su4)./(4*cosd(theta)*sind(theta));
                    
                    %need to compute omega, the radian frequency
                    omega = 2*pi.*fu;dw = omega(3)-omega(2);
                    
                    %by definition, the covariance of two signals is the integral of
                    %the cospectrum (this bit is from Gerbi, 2008 eqn. 6)
                    uwstar = trapz(COuw);vwstar = trapz(COvw);
                    
                    %convert frequencies to wavenumbers via the Taylor Hypothesis
                    kcu = 2*pi*cf./abs(nanmean(X));kcv = 2*pi*cf./abs(nanmean(Y)); %cutoff k
                    ku = omega./abs(nanmean(X));kv = omega./abs(nanmean(Y));           %wavenumber (k)
                    
                    %compute variance-preserving cospectrum
                    COuwvar = COuw.*ku;
                    COvwvar = COvw.*kv;
                    
                    %calculate model cospectrum
                    pk = max(findpeaks(COuwvar));
                    ku0 = ku(COuwvar == pk);
                    pk = max(findpeaks(COvwvar));
                    kv0 = kv(COvwvar == pk);
                    COuwstar = uwstar*((7/3*pi)*sin(3*pi/7))*((1/ku0)./(1+(ku./ku0).^(7/3)));
                    COvwstar = vwstar*((7/3*pi)*sin(3*pi/7))*((1/kv0)./(1+(ku./kv0).^(7/3)));
                    
                    %compute Ogive curve (CDF) of cospectra
                    oguw = cumtrapz(COuw);
                    ogvw = cumtrapz(COvw);
                    oguwstar = cumtrapz(COuwstar.*dw);
                    ogvwstar = cumtrapz(COvwstar.*dw);
                    
                    %turb energy is contained in low f in CDF functions
                    %lin fit the ogive curves, intercept is RS
                    %%%CHANGE THIS LINE TO USE THE \ OPERATOR!! 
                    pf = polyfit(oguwstar(ku < kcu),oguw(ku < kcu),1);
                    uw = pf(2);
                    pf = polyfit(ogvwstar(kv < kcv),ogvw(kv < kcv),1);
                    vw = pf(2);
                    
                    if sign(uwstar) ~= sign(uw);
                        uw = uw*-1;
                    end
                    if sign(vwstar) ~= sign(vw);
                        vw = vw*-1;
                    end
                    
                    %save variables in structures
                    %reject instances where rolloff k is greater than cutoff k
                    if ku0 > round(kcu)
                        UW(jj) = NaN;
                    else
                        UW(jj) = uw;
                    end
                    if kv0 > round(kcv)
                        VW(jj) = NaN;
                    else
                        VW(jj) = vw;
                    end
                    UWSTAR(jj) = uwstar;
                    VWSTAR(jj) = vwstar;
                    KU0(jj) = ku0;
                    KV0(jj) = ku0;
                end
            end
            RS.(date).(fn{ii}).time(j,:) = time2;
            RS.(date).(fn{ii}).uw(j,:) = UW;
            RS.(date).(fn{ii}).vw(j,:) = VW;
            RS.(date).(fn{ii}).uwstar(j,:) = UWSTAR;
            RS.(date).(fn{ii}).vwstar(j,:) = VWSTAR;
            RS.(date).(fn{ii}).ku0(j,:) = KU0; %for estimating horizontal length scales (see: Kirincich 2010).
            RS.(date).(fn{ii}).kv0(j,:) = KV0;
        end
        matlabpool CLOSE
    end
end
disp(['Processing time ' num2str(toc/60) ' minutes'])
savedatdir = 'd:\Projects\Mekong_W2015\DataAnalysis\Paper2\';
fname = 'RStress_10min';
if orienttest
    fname = [fname '+_-5deg'];
end
save([savedatdir fname],'RS','-v7.3')