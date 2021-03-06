%Calculate Reynolds stress by calculating the cospectra of velocities, then
%filtering by the turbulent frequency band previously defined by our other
%analyses. This is version 5 of this script.
clear
close all
datdir = 'd:\Projects\Mekong_W2015\DataAnalysis\Paper3\VPs\';
files = {'7March2015_Vels.mat';'8March2015_Vels.mat';'10March2015_Vels.mat';...
    '14March2015a_Vels.mat'};
RS = struct();
days = {'day1';'day2';'day3';'day4'};
%ENABLE THESE IF YOU WANT TO MANUALLY ROTATE THE VELOCITIES
orienttest = 1;
angle = -0.5; %rotation angle for testing limitations
pitch = 0;
roll = 1;
%These are from PCA of the velocities to define the mean current direction
heading = [55 240 55 260]; %headings are in order! (HTA1 - 3, VTA) was 55
cutoff = [1.4 1.35 1.04 1.09]; %from wave cutoff freq (see CutoffFreqs)
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
        cf = cutoff(i);
        disp(['Applying cutoff: ' sprintf('%0.2f',cf) 'Hz'])
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
        disp(['Loading & rotations completed at: ' datestr(now,'HH:MM:SS')])
        if orienttest
            %load in k0 data from original run if running orientation test
            kzed = load('D:\Projects\Mekong_W2015\DataAnalysis\Paper2\RS_estimates\RS_k0_values.mat');
            %perform an additional rotation to test if RS changes by
            %altering the angles by a small amount
            disp('User requests +/-0.5 deg pitch or roll adjustment to velocities')
            if pitch
                pp = angle*pi/180;
                rr = 0;
            elseif roll
                pp = 0;
                rr = angle*pi/180;
            end
            P = [cos(pp) -sin(pp)*sin(rr) -cos(rr)*sin(pp);...
                0             cos(rr)          -sin(rr);  ...
                sin(pp) sin(rr)*cos(pp)  cos(pp)*cos(rr)];
            x = zeros(m,n);y = zeros(m,n);z1 = zeros(m,n);z2 = zeros(m,n);
            for k = 1:n
                V = [b1(:,k) b2(:,k) b3(:,k)];
                vels = V*P;
                x(:,k) = vels(:,1);
                y(:,k) = vels(:,2);
                z1(:,k) = vels(:,3);
                z2(:,k) = z2(:,k); %4-beam ADCP transforms ignore beam 4.
            end
            fprintf('Applied %0.1f deg pitch adjustment and %0.1f roll adjustment to velocities\n',pp*180/pi,rr*180/pi)
        else
            x = b1;y = b2;z1 = b3;z2 = b4;
        end
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
            %Define variables to save out from jj loop
            UW = zeros(1,35);
            VW = zeros(1,35);
            UWSTAR = zeros(1,35);
            VWSTAR = zeros(1,35);
            UWINT = zeros(1,35);
            VWINT = zeros(1,35);
            KU0 = zeros(1,35);
            KV0 = zeros(1,35);
            FU0 = zeros(1,35);
            FV0 = zeros(1,35);
            ubrU = zeros(1,35);
            ursq = zeros(1,35);
            vrsq = zeros(1,35);
            K = zeros(4097,35);
            KC = zeros(1,35);
            cspecuw = zeros(4097,35);
            cspecvw = zeros(4097,35);
            cspecuws = zeros(4097,35);
            cspecvws = zeros(4097,35);
            oguw = zeros(4097,35);
            ogvw = zeros(4097,35);
            oguws = zeros(4097,35);
            ogvws = zeros(4097,35);
            time2 = time(ind(j));
            for jj = 1:35
                X = x(idx,jj);
                Y = y(idx,jj);
                Z1 = z1(idx,jj);
                Z2 = z2(idx,jj);
                if min(X) == 0 && max(X) == 0 || min(Y) == 0 && max(Y) == 0 %skip bins with zero velocity
                    UW(jj) = NaN;
                    VW(jj) = NaN;
                    UWSTAR(jj) = NaN;
                    VWSTAR(jj) = NaN;
                    UWINT(jj) = NaN;
                    VWINT(jj) = NaN;
                    KU0(jj) = NaN;
                    KV0(jj) = NaN;
                    K(:,jj) = NaN(4097,1);
                    KC(jj) = NaN;
                    cspecuw(:,jj) = NaN(4097,1); %observation
                    cspecvw(:,jj) = NaN(4097,1);
                    cspecuws(:,jj) = NaN(4097,1); %model
                    cspecvws(:,jj) = NaN(4097,1);
                    oguw(:,jj) = NaN(4097,1); %observation
                    ogvw(:,jj) = NaN(4097,1);
                    oguws(:,jj) = NaN(4097,1); %model
                    ogvws(:,jj) = NaN(4097,1);
                else
                    %Update: 28/03/18 Steve recommends rotating vertical velocities to
                    %the principal component direction for the near-bed experiments
                    %for future reference, plug the whole t-s into cmgpca
                    %and use a single rotation value for all of the data.
                    if i == 1 || i == 4
                        %                         out = cmgpca(X,Z1);
                        thz1 = ((2))*pi/180;
                        Rz1 = [cos(thz1) -sin(thz1); sin(thz1) cos(thz1)];
                        %                         out = cmgpca(X,Z2);
                        thz2 = ((2))*pi/180;
                        Rz2 = [cos(thz2) -sin(thz2); sin(thz2) cos(thz2)];
                        xz1 = [X Z1]*Rz1;
                        xz2 = [X Z2]*Rz2;
                        X = xz1(:,1);
                        Z1 = xz1(:,2);
                        Z2 = xz2(:,2);
                    end
                    u = detrend(X);
                    v = detrend(Y);
                    w1 = detrend(Z1);
                    w2 = detrend(Z2);
                    %estimate ratio of wave orbital velocities to mean
                    %current velocities
                    [Guu,f] = pwelch(u,[],[],[],fs);
                    [Gvv,~] = pwelch(v,[],[],[],fs);
                    Guv = Guu + Gvv;
                    range = find(f>0.05&f<cf);
                    df = f(3)-f(2);
                    ubr = sqrt(2*sum(Guv(range)*df)); %wave orbital velocity in waveband
                    U = mean(sqrt((u.^2)+(v.^2))); %avg. horizontal speed
                    uwU = ubr/U; %ratio of orbital velocity to mean current speed
                    ubrU(jj) = uwU;
                    if uwU > 2 %reject bins where this ratio is large
                        UW(jj) = NaN;
                        VW(jj) = NaN;
                        UWSTAR(jj) = NaN;
                        VWSTAR(jj) = NaN;
                        UWINT(jj) = NaN;
                        VWINT(jj) = NaN;
                        KU0(jj) = NaN;
                        KV0(jj) = NaN;
                        K(:,jj) = NaN(4097,1);
                        cspecuw(:,jj) = NaN(4097,1); %observation
                        cspecvw(:,jj) = NaN(4097,1);
                        cspecuws(:,jj) = NaN(4097,1); %model
                        cspecvws(:,jj) = NaN(4097,1);
                        oguw(:,jj) = NaN(4097,1); %observation
                        ogvw(:,jj) = NaN(4097,1);
                        oguws(:,jj) = NaN(4097,1); %model
                        ogvws(:,jj) = NaN(4097,1);
                    else
                        theta = 30;
                        %calc along-beam velocities
                        u1 = -u*sind(theta)-w1*cosd(theta);
                        u2 = u*sind(theta)-w2*cosd(theta);
                        u3 = -v*sind(theta)-w1*cosd(theta);
                        u4 = v*sind(theta)-w2*cosd(theta);
                        
                        %compute velocity spectra
                        [Su1,fu,pu1] = pwelch(u1,[],[],[],fs,'confidencelevel',0.95);
                        [Su2,~,pu2] = pwelch(u2,[],[],[],fs,'confidencelevel',0.95);
                        [Su3,~] = pwelch(u3,[],[],[],fs);
                        [Su4,~] = pwelch(u4,[],[],[],fs);
                        
                        %calculate cospectra
                        COuw = (Su1-Su2)./(4*cosd(theta)*sind(theta));
                        COvw = (Su3-Su4)./(4*cosd(theta)*sind(theta));
                        
                        %Edit: 22/05/18
                        %flip cospectra u/d so to focus only on the hf
                        %turbulence now BELOW the waveband (after flipping)
                        COuw = flipud(COuw);
                        COvw = flipud(COvw);
                        
                        %cospectral estimates of error bounds
                        COuwerrl = (pu1(:,1)-pu2(:,1))./(4*cosd(theta)*sind(theta));
                        COuwerrh = (pu1(:,2)-pu2(:,2))./(4*cosd(theta)*sind(theta));
                        COuwerrl = flipud(COuwerrl);
                        COuwerrh = flipud(COuwerrh);
                        
                        %estimate wavenumbers based on U
                        k = 2*pi*fu./U;kc = 2*pi*cf/U;
                        k = flipud(k); %also flip k (see edit: 22/05/18)
                        dk = k(2)-k(3);
                        df = 2*pi*fu(3)-2*pi*fu(2);
                        
                        %by definition, the covariance of two signals is the integral of
                        %the cospectrum (this bit is from Gerbi, 2008 eqn. 6)
                        uwstar = sum(COuw.*df);
                        vwstar = sum(COvw.*df);
                        
                        %compute variance-preserving cospectrum
                        COuwvar = COuw.*k;
                        COvwvar = COvw.*k;
                        
                        %if running the rotations for QC, load in the ku0
                        %and kv0 values generated in the first run. DO NOT
                        %USE THIS OPTION IF YOU ARE RUNNING THE FIRST RUN.
                        if orienttest
                            %compute Ogive curve (CDF) of cospectra
                            oguw = cumsum(COuw.*dk);
                            ogvw = cumsum(COvw.*dk);
                            %use these k's as "default" values
                            ku0 = kzed.(date).(fn{ii}).ku0(j,jj);
                            kv0 = kzed.(date).(fn{ii}).kv0(j,jj);
                            kc = kzed.(date).(fn{ii}).kc(j,jj);
                            A = (7/(3*pi))*sin(3*pi/7);
                            B = (1/ku0)./(1+((k./ku0).^(7/3)));
                            COuwstar = (A*B);
                            B = (1/kv0)./(1+((k./kv0).^(7/3)));
                            COvwstar = (A*B);
                            %compute model ogives
                            oguwstar = cumsum(COuwstar.*dk);
                            ogvwstar = cumsum(COvwstar.*dk);
                            %calculate stresses
                            [uw,~,rsq] = my_regression(oguw(k>kc),oguwstar(k>kc),1);
                            ursq(jj) = rsq;
                            uw = uw(1);
                            
                            [vw,~,rsq] = my_regression(ogvw(k>kc),ogvwstar(k>kc),1);
                            vrsq(jj) = rsq;
                            vw = vw(1);
                        else
                            %To estimate k0, the dominant length scale of the
                            %stress-carrying eddies, we need to iteratively
                            %solve for it using uw and the r-squared from the
                            %ogive fits. Steve recommends solving for uw and
                            %r-squared over the range of k >= kc.
                            ktest = linspace(kc,max(k),40);
                            ursqt = zeros(1,40);
                            vrsqt = zeros(1,40);
                            uwt = zeros(1,40);
                            vwt = zeros(1,40);
                            %compute Ogive curve (CDF) of cospectra
                            oguw = cumsum(COuw.*dk);
                            ogvw = cumsum(COvw.*dk);
                            for qq = 1:40
                                k0 = ktest(qq);
                                A = (7/(3*pi))*sin(3*pi/7);
                                B = (1/k0)./(1+((k./k0).^(7/3)));
                                COuwstar = (A*B);
                                COvwstar = (A*B);
                                %compute the model with each successive k0
                                oguwstar = cumsum(COuwstar.*dk);
                                ogvwstar = cumsum(COvwstar.*dk);
                                
                                %turb energy is contained in low f in CDF functions
                                %lin fit the ogive curves, Edit: 23/05/2018 Steve
                                %recommends fitting through zero; instead of using
                                %polyfit (which provides an offset, the y-int), use
                                %"my_regression" which fits through zero.
                                
                                [uw,~,rsq] = my_regression(oguw(k>kc),oguwstar(k>kc),1);
                                ursqt(qq) = rsq;
                                uwt(qq) = uw(1);
                                
                                [vw,~,rsq] = my_regression(ogvw(k>kc),ogvwstar(k>kc),1);
                                vrsqt(qq) = rsq;
                                vwt(qq) = vw(1);
                            end
                            %now estimate max rsq for uw and vw. Save k0 for
                            %each and filter the results
                            [~,rsid] = max(ursqt);
                            ku0 = ktest(rsid);
                            uw = uwt(rsid);
                            ursq(jj) = ursqt(rsid);
                            
                            [~,rsid] = max(vrsqt);
                            kv0 = ktest(rsid);
                            vw = vwt(rsid);
                            vrsq(jj) = vrsqt(rsid);
                            
                            %recompute model, ogives
                            A = (7/(3*pi))*sin(3*pi/7);
                            B = (1/ku0)./(1+((k./ku0).^(7/3)));
                            COuwstar = (A*B);
                            B = (1/kv0)./(1+((k./kv0).^(7/3)));
                            COvwstar = (A*B);
                            oguwstar = cumsum(COuwstar.*dk);
                            ogvwstar = cumsum(COvwstar.*dk);
                        end
                        %Compute above-waveband integrated uw and vw (Gerbi
                        %et al. 2008). Should be similar, but smaller than
                        %uw and vw computed by fitting the model
                        uwint = sum(COuw(k > kc).*dk);
                        vwint = sum(COvw(k > kc).*dk);
                        
                        if sign(uwstar) ~= sign(uw);
                            uw = uw*-1;
                        end
                        if sign(vwstar) ~= sign(vw);
                            vw = vw*-1;
                        end
                        if sign(uwstar) ~= sign(uwint);
                            uwint = uwint*-1;
                        end
                        if sign(vwstar) ~= sign(vwint);
                            vwint = vwint*-1;
                        end
                        
                        %save variables in structures
                        %reject instances where rolloff k is less than
                        %cutoff k, or where r-squared values are less than
                        %0.7
                        if ku0 <= kc || ursq(jj) < 0.7 %see Davis & Monismith (2011)
                            UW(jj) = NaN;
                            UWINT(jj) = NaN;
                            cspecuw(:,jj) = NaN(4097,1);
                            cspecuws(:,jj) = NaN(4097,1);
                            oguw(:,jj) = NaN(4097,1);
                            oguws(:,jj) = NaN(4097,1);
                        else
                            UW(jj) = uw;
                            UWINT(jj) = uwint;
                            cspecuw(:,jj) = COuw;
                            cspecuws(:,jj) = COuwstar;
                            oguw(:,jj) = oguw;
                            oguws(:,jj) = oguwstar;
                        end
                        if kv0 <= kc || vrsq(jj) < 0.7
                            VW(jj) = NaN;
                            VWINT(jj) = NaN;
                            cspecvw(:,jj) = NaN(4097,1);
                            cspecvws(:,jj) = NaN(4097,1);
                            ogvw(:,jj) = NaN(4097,1);
                            ogvws(:,jj) = NaN(4097,1);
                        else
                            VW(jj) = vw;
                            VWINT(jj) = vwint;
                            cspecvw(:,jj) = COvw;
                            cspecvws(:,jj) = COvwstar;
                            ogvw(:,jj) = ogvw;
                            ogvws(:,jj) = ogvwstar;
                        end
                        UWSTAR(jj) = uwstar;
                        VWSTAR(jj) = vwstar;
                        KU0(jj) = ku0;
                        KV0(jj) = kv0;
                        FU0(jj) = ku0*U/2*pi;
                        FV0(jj) = kv0*U/2*pi;
                        K(:,jj) = k;
                        KC(jj) = kc;
                    end
                end
            end
            RS.(date).(fn{ii}).time(j,:) = time2;
            RS.(date).(fn{ii}).uw(j,:) = UW;
            RS.(date).(fn{ii}).vw(j,:) = VW;
            RS.(date).(fn{ii}).uwstar(j,:) = UWSTAR;
            RS.(date).(fn{ii}).vwstar(j,:) = VWSTAR;
            RS.(date).(fn{ii}).uwint(j,:) = UWINT;
            RS.(date).(fn{ii}).vwint(j,:) = VWINT;
            RS.(date).(fn{ii}).ku0(j,:) = KU0; %for estimating horizontal length scales (see: Kirincich 2010).
            RS.(date).(fn{ii}).kv0(j,:) = KV0;
            RS.(date).(fn{ii}).fu0(j,:) = FU0;
            RS.(date).(fn{ii}).fv0(j,:) = FV0;
            RS.(date).(fn{ii}).k(j,:) = mean(K,2);
            RS.(date).(fn{ii}).kc(j,:) = KC;
            RS.(date).(fn{ii}).COuw(j,:) = nanmean(cspecuw,2);
            RS.(date).(fn{ii}).COvw(j,:) = nanmean(cspecvw,2);
            RS.(date).(fn{ii}).COuwstar(j,:) = nanmean(cspecuws,2);
            RS.(date).(fn{ii}).COvwstar(j,:) = nanmean(cspecvws,2);
            RS.(date).(fn{ii}).oguw(j,:) = nanmean(oguw,2);
            RS.(date).(fn{ii}).ogvw(j,:) = nanmean(ogvw,2);
            RS.(date).(fn{ii}).oguwstar(j,:) = nanmean(oguws,2);
            RS.(date).(fn{ii}).ogvwstar(j,:) = nanmean(oguws,2);
            RS.(date).(fn{ii}).ubrU(j,:) = ubrU;
            RS.(date).(fn{ii}).ursq(j,:) = ursq;
            RS.(date).(fn{ii}).vrsq(j,:) = vrsq;
        end
        clear dat
    end
    disp(['Processing time ' num2str(toc/60) ' minutes'])
end
disp(['Files completed at: ' datestr(now)])
fprintf('End of Run\n\n')
savedatdir = 'd:\Projects\Mekong_W2015\DataAnalysis\Paper2\RS_estimates\';
fname = 'RStress_10min_v5';
if orienttest
    if pitch
        txt = 'pitch';
    elseif roll
        txt = 'roll';
    end
    fname = [fname '_' num2str(sign(angle)) 'deg_' txt];
end
save([savedatdir fname],'RS','-v7.3')