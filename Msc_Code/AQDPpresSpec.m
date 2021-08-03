%Plot wave spectra for Karin from the F2F Deployment. Plot Spectra from
%(most of) the pressure sensors along the transect line. This script will
%be only for plotting Aquadopps, since they have a different sample rate
%from the pressure loggers.

clear
filelist = {'HR3_7March2015_f_pad.mat';'AD5116_9March2015_f_pad.mat';...
    'AD5117_9March2015_f_pad.mat'};
name = {'HR3';'AD5116';'AD5117'};
fdir = 'C:\Users\bkn5\Projects\Mekong_W2015\Data\Aquadopp\F2F2\';
dname = 'F2F2'; %deployment
plotfigs = 0;
savefigs = 0;
z = 0.075; %pressure sensor height
intv = 10; %averaging interval (minutes)
start = datenum(2015,03,06,14,24,00); %pick a time
stop = datenum(2015,03,06,14,34,00);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for kk = 1:length(filelist)
    disp(['Loading file ' filelist{kk}])
    load([fdir filelist{kk}])
    %since the AQDPs recorded over multiple cycles, provide limits for when the
    %instrument was underwater.
    ind = find(aqdp.datenum >= start & aqdp.datenum <=stop);
    datetime = aqdp.datenum(ind);
    pressure = aqdp.pressure(ind);
    temperature = aqdp.temperature(ind);
    
    %Initialize Global Variables
    fs = str2double(aqdp.metadata.samprate(1:2)); %sampling frequency
    avt = fs*intv*60; %# samples in interval
    idx = avt:avt:length(pressure);
    idxx = [1 idx];
    it = 1;
    while it < length(idxx)
        disp(['Iteration number ' num2str(it)])
        ind = idxx(it):idxx(it+1);
        p = pressure(ind);
        t = temperature(ind);
        if rem(length(p), 2) == 1
            p = p(1:end-1); %force inputs to be even
            t = t(1:end-1);
        end
        %fill gaps in temperature/pressure signal
        if cmgidgaps(p) > 0
            disp(['Found ' num2str(cmgidgaps(p)) ' gaps in pressure time-series'])
            disp(['Found ' num2str(cmgidgaps(t)) ' gaps in temperature time-series'])
            nlin = 1E2;
            maxgaps = 1E5;
            p = cmgbridge(p,nlin,maxgaps,maxgaps);
            t = cmgbridge(t,nlin,maxgaps,maxgaps);
            if cmgidgaps(p) > 0 || cmgidgaps(t) > 0
                p(isnan(p)) = 0;
                t(isnan(p)) = 0;
                disp(['Number of gaps in pressure remaining: ' num2str(cmgidgaps(p))])
                disp(['Number of gaps in temperature remaining: ' num2str(cmgidgaps(t))])
            else
                disp('Gaps filled')
            end
        end
        salt = repmat(20,length(t),1); %assume constant salinity of 20PSU
        rho = SeaDensity(salt,t,p);
        
        %Calculate Depth (h) from the pressure signal %%%%%%%%%%%%%%%%%%%%%
        %designate start/stop time
        %From: UNESCO (1983): Algorithms for computation of fundamental properties
        %of seawater. UNESCO technical papers in marine science 44:1-55.
        x = (sin(aqdp.metadata.lat/57.29578))^2;
        g = zeros(length(p),1);h = zeros(length(p),1);
        for i = 1:length(p)
            g(i,:) = 9.780318*(1.0+(5.2788E-3+2.36E-5*x)*x)+1.092E-6*p(i,:);
            h(i,:) = ((((-1.82E-15*p(i,:)+2.279E-10)*p(i,:)-2.2512E-5)*p(i,:)+9.72659)*p(i,:))/g(i,:);
        end
        h = mean(h);
        g = mean(g);
        rho = mean(rho);
        
        %%2a. Filter timeseries to remove low-f and higher f-noise %%%%%%%%%%%%
        nf = fs/2; %nyquist criterium
        %%2b. Compute fft for spectral analysis of pressure data %%%%%%%%%%%%%%
        ps = detrend(p);
        n = length(ps);
        pf = fft(ps,length(p));
        freq = 0:fs/length(ps):fs/2;
        
        pf = pf(1:n/2+1); %get fft(ps) corresponding to [0 nf)
        if plotfigs
            % plot results
            figure
            plot(freq(1:length(p)/2+1),2*abs(pf(1:length(p)/2+1)))
            xlabel('Frequency (Hz)')
            ylabel('Amplitude (fft(p))')
            set(gca,'xlim',[0 fs])
            prompt = 'Is the filtering acceptable? [y/n] ';
            result = input(prompt,'s');
            if strcmp(result,'y') || strcmp(result,'yes');
            end
            if strcmp(result,'n') || strcmp(result,'no');
                return
            end
        end
        
        %%3a. Calculate wavenumber (k) %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        disp('Computing surface pressure attenuation adjustment')
        k = zeros(length(freq),1);
        attn = zeros(length(freq),1); %attenuation
        T = 1./freq; %period = 1/freq
        omega = 2*pi./T';
        for i = 1:length(freq)
            k(i) = om2k(omega(i),h);
            %%3b. Calculate pressure attentuation for each f %%%%%%%%%%%%%%%%%%
            %From Raubenheimer, B et al. 1998 "Estimating Wave Heights from
            %Pressure Measured in Sand Bed"
            attn(i) = (cosh(k(i)*(z)))/(cosh(k(i)*h));
        end
        %%3c. Calculate attenuated surface pressure %%%%%%%%%%%%%%%%%%%%%%%%%%%
        %for the ifft to be real, it must have Hermetian Symmetry.
        disp('Applying adjustment to timeseries')
        pfa = pf./attn; %adjusted fourier transfm.
        pfa(attn<0.33) = 0; %zero out attenuations below 0.33Hz

%         %matlab's fft functions are not normalized.
%         pf = pf/(sqrt(numel(pf)));
%         %check to see this meets parseval's theorem
%         if roundn(norm(ps),-8) == roundn(norm(pf),-8)
%             disp('Normalization of fft(pressure) meets Parsevals Theorem')
%         else
%             disp('WARNING: normalization of fft(pressure) was unsuccessful')
%             disp(['t-s: ' sprintf('%.10f',norm(ps))])
%             disp(['fft(t-s): ' sprintf('%.10f',norm(pf))])
%             return 
%         end
        
        psurf = ifft(pfa,length(ps),'symmetric');
        psdpf = (1/(fs*n))*abs(pfa).^2;
        psdpf(2:end-1) = 2*psdpf(2:end-1); %psd of ps
        Spec.(name{kk}).datetime(it,:) = datetime(1:end-1);
        Spec.(name{kk}).psurf(it,:) = psurf;
        Spec.(name{kk}).p(it,:) = p;
        Spec.(name{kk}).ps(it,:) = ps;
        Spec.(name{kk}).psd(it,:) = psdpf;
        Spec.(name{kk}).freq = freq;
        it = it+1;
    end
    clear aqdp
end
disp('Finished!')
save AqdpSpectrum.mat Spec