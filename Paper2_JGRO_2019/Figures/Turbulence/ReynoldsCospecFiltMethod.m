%Calculate Reynolds stress by calculating the cospectra of velocities, then
%filtering by the turbulent frequency band previously defined by our other
%analyses. 

clear
datdir = 'c:\Users\bkn5\Projects\Mekong_W2015\DataAnalysis\Paper1\';
load([datdir 'HTAday1Vels.mat'])
ploton = 0;
cutoff = 6;                                                         %Hz, this is the cutoff freq from the spectral slopes
x = HTA.vpro1.x;y = HTA.vpro1.y;z1 = HTA.vpro1.z1;z2 = HTA.vpro1.z2;
%subdivide timeseries into sections 10 minutes long
intv = 10;
fs = 50;
avt = intv*fs*60;                                                   %samples/interval
ind = [1 avt:avt:length(x)];
for i = 1:length(ind)-1                                             %loop in time
    for ii = 1:25                                                   %loop by depth bin
        X = x(ind(i):ind(i+1),ii);
        Y = y(ind(i):ind(i+1),ii);
        Z1 = z1(ind(i):ind(i+1),ii);
        Z2 = z2(ind(i):ind(i+1),ii);
        %force inputs to be even
        if rem(length(X),2) == 1
            X = X(1:end-1,:);
        end
        if rem(length(Y),2) == 1
            Y = Y(1:end-1,:);
        end
        if rem(length(Z1),2) == 1
            Z1 = Z1(1:end-1,:);
        end
        if rem(length(Z2),2) == 1
            Z2 = Z2(1:end-1,:);
        end
        %compute along-beam velocities
        u = detrend(X);
        v = detrend(Y);
        w1 = detrend(Z1);
        w2 = detrend(Z2);
        theta = 30;                                                  %deg; half angle of beams to the vertical
        
        u1 = -u*sind(theta)-w1*cosd(theta);
        u2 = u*sind(theta)-w2*cosd(theta);
        u3 = -v*sind(theta)-w1*cosd(theta);
        u4 = v*sind(theta)-w2*cosd(theta);
        
        %compute velocity spectra
        [m,~] = size(u);
        nfft = 0.25*m;
        window = hanning(m,'periodic');
        noverlap = 0.7*length(window);
        
        [Su1,fu] = pwelch(u1,window,noverlap,nfft,fs);
        [Su2,~] = pwelch(u2,window,noverlap,nfft,fs);
        [Su3,~] = pwelch(u3,window,noverlap,nfft,fs);
        [Su4,~] = pwelch(u4,window,noverlap,nfft,fs);
        
        %calculate cospectra
        COuw = (Su1-Su2)./(4*cosd(theta)*sind(theta));
        COvw = (Su3-Su4)./(4*cosd(theta)*sind(theta));
    
        %integrate cospectra: this is the estimate of Reynolds stress
        %before filtering (the covariance)
        uwstar = trapz(COuw);vwstar = trapz(COvw);
        
        %convert frequencies to wavenumbers via the Taylor Hypothesis
        kcu = cutoff./abs(nanmean(X));kcv = cutoff./abs(nanmean(Y)); %cutoff k
        ku = fu./abs(nanmean(X));kv = fu./abs(nanmean(Y));           %wavenumber (k)
        
        %Generate model spectra:
        %Compute k0 from variance-preserving cospectra
        cospectu = COuw.*ku;                                         %variance-preserving cospectrum
        fun = fit(ku,cospectu,'gauss1',...                           %compute gaussian fit to represent the cospectrum
            'normalize','on');
        yucoef = coeffvalues(fun);
        yu = fun(ku);
        [~,id] = findpeaks(abs(yu));
        k0uw = ku(id);                                               %k0 is the peak of the variance-preserving cospectrum (observations)
        
        cospectv = COvw.*kv;
        fun = fit(kv,cospectv,'gauss1',...                           
            'normalize','on');
        yvcoef = coeffvalues(fun);
        yv = fun(kv);
        [~,id] = findpeaks(abs(yv));                                     
        k0vw = kv(id);

        if ploton                                                    %(optional) plot var-preserving cospectra gaussian model fits;
            figure                                                   %peak of the model is the observed estimate for k0, the rolloff wn
            subplot(211)
            semilogx(ku,cospectu,'o'),hold on
            semilogx(ku,yu,'--r','LineWidth',1.5)
            xlabel('Wavenumber (k)'), ylabel('kCo_u_w*(k) (rad ms^-^2)')
            subplot(212)
            semilogx(kv,cospectv,'o'),hold on
            semilogx(kv,yv,'--r','LineWidth',1.5)
            xlabel('Wavenumber (k)'), ylabel('kCo_v_w*(k) (rad ms^-^2)')
            suptitle(['Burst ' num2str(i) ', Bin #' num2str(ii)])
        end
        
        if ~(k0uw > kcu)                                             %cutoff criteria; k0 must be greater than the cutoff k
            
            %Generate model cospectra from covariance and k0
            COuwstar = -uwstar*((7/3*pi)*sin(3*pi/7))*((1/k0uw)./(1+(ku./k0uw).^(7/3)));
            [~,id] = findpeaks(abs(COuwstar.*ku));
            k0uwstar = ku(id);
            
            %compute Ogive curves
            OGuw = cumtrapz(COuw);
            OGuwstar = cumtrapz(COuwstar);
            
            %calculate Reynolds stress
            id = find(ku > kcu);                                     %find wavenumbers below the cutoff wn
            RSuw = trapz(COuw(id));
            RSuwstar = trapz(COuwstar(id));
            
            if ploton
                figure
                subplot(311)
                semilogx(ku./k0uw,COuwstar,'LineWidth',2)
                set(gca,'XTickLabel',[],'Xlim',[1E-2 1E2])
                ylabel('Co_u_w*(k) (m^2s^-^2)')
                subplot(312)
                semilogx(ku./k0uw,COuwstar.*ku,'LineWidth',2)
                set(gca,'XTickLabel',[],'Xlim',[1E-2 1E2])
                ylabel('kCo_u_w*(k) (rad ms^-^2)')
                subplot(313)
                semilogx(ku./k0uw,OGuwstar,'LineWidth',2)
                set(gca,'Xlim',[1E-2 1E2])
                ylabel('\int Co_u_w*(k) (m^2s^-^2)')
                xlabel('k/k_0')
                suptitle(['Burst ' num2str(i) ', Bin #' num2str(ii)])
                
                figure
                semilogx(ku(id),OGuw(id),'-k','LineWidth',1.5), hold on
                semilogx(ku(id),OGuwstar(id),'--','Color',[0.5 0.5 0.5],'LineWidth',1.5)
                ylabel('\int Co_u_wdk (m^2s^-^2)')
                xlabel('k (rad/m)')
            end
            
            %save variables in structure
            rss.obs.uw(i,ii) = RSuw;
            rss.obs.k0uw(i,ii) = k0uw;
            rss.mod.uw(i,ii) = RSuwstar;
            rss.mod.k0uw(i,ii) = k0uwstar;
            
        elseif ~(k0uw < kcu)
            %save variables in structure
            rss.obs.uw(i,ii) = 0;
            rss.obs.k0uw(i,ii) = 0;
            rss.mod.uw(i,ii) = 0;
            rss.mod.k0uw(i,ii) = 0;
        end
       if ~(k0vw > kcv)                                              %cutoff criteria; k0 must be greater than the cutoff k
            
            %Generate model cospectra from covariance and k0
            COvwstar = vwstar*((7/3*pi)*sin(3*pi/7))*((1/k0vw)./(1+(kv./k0vw).^(7/3)));
            [~,id] = findpeaks(abs(COvwstar.*kv));
            k0vwstar = kv(id);
            
            %calculate Reynolds stress
            id = find(kv > kcv);                                     %find wavenumbers below the cutoff wn
            RSvw = trapz(COvw(id));
            RSvwstar = trapz(COvwstar(id));
            
            %save variables in structure
            rss.obs.vw(i,ii) = RSvw;
            rss.obs.k0vw(i,ii) = k0vw;
            rss.mod.vw(i,ii) = RSvwstar;
            rss.mod.k0vw(i,ii) = k0vwstar;
            
        elseif ~(k0vw < kcv)
            %save variables in structure
            rss.obs.vw(i,ii) = 0;
            rss.obs.k0vw(i,ii) = 0;
            rss.mod.vw(i,ii) = 0;
            rss.mod.k0vw(i,ii) = 0;
       end
    end
end

        
        
