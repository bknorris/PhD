%This script calculates Wave Power Flux (WPF) from Aquadopp pressure data,
%and differences WPF with TKE dissipation measurements collected by our
%Vectrinos during the same experiment over different days.

%%Preliminary Checklist:
%Terms (units):
%k = wave number
%z = instrument hab (m) [nominal height of pressure sensor]
%h = water depth (m)
%H = sig. wave height (m)
%omega = wave angular velocity (rad*m/s)
%g = gravity (m/s^2)
%rho = water density (kg/m^3)
%Cg = wave group velocity (m/s)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% clear %a good idea
load('c:\Users\bkn5\Projects\Mekong_W2015\Data\Aquadopp\FSS\HR3_9March2015_f_pad.mat')
dname = 'FSS32'; %deployment
name = 'ADHR3'; %instrument
plotfigs = 0;
plotwpf = 0;
z = 0.075; %pressure sensor height
intv = 10; %averaging interval (minutes)
start = datenum(2015,03,8,14,14,38);
stop = datenum(2015,03,8,16,40,38);
savedatfile = 1;
% plot(aqdp.datenum,aqdp.pressure)
% datetickzoom('x','keepticks','keeplimits')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%since the AQDPs recorded over multiple cycles, provide limits for when the
%instrument was underwater.
ind = find(aqdp.datenum >= start & aqdp.datenum <=stop);
pressure = aqdp.pressure(ind);
temperature = aqdp.temperature(ind);
%Initialize Global Variables %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
WPF = struct();
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
    
    %%1. Calculate Depth (h) from the pressure signal %%%%%%%%%%%%%%%%%%%%%
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
%     [b,a] = butter(5,0.08/nf,'high');
%     p = filtfilt(b,a,p); %filter pressure data
%     [b,a] = butter(5,2/nf,'low');
%     p = filtfilt(b,a,p); %filter pressure data
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
        attn(i) = (cosh(k(i)*(h+z)))/(cosh(k(i)*h));
    end
    if z > 0
        attn(attn < 0.2) = 0.2; %eliminate over-reduction less than 0.02Hz
        attn(freq<0.01 | freq>(nf/1.5)) = 1; % only for selected range (0.01-2.67Hz)
    end
    %%3c. Calculate attenuated surface pressure %%%%%%%%%%%%%%%%%%%%%%%%%%%
    %for the ifft to be real, it must have Hermetian Symmetry.
    disp('Applying adjustment to timeseries')
    pfa = pf.*attn; %adjusted fourier transfm.    
    psurf = ifft(pfa,length(ps),'symmetric');
%     psurf = smooth(psurf,10,'loess');

    %matlab's fft functions are not normalized.
    pfa = pfa/(sqrt(numel(pfa)));
    %check to see this meets parseval's theorem
    if roundn(norm(ps),-8) == roundn(norm(pfa),-8)
        disp('Normalization of fft(pressure) meets Parsevals Theorem')
    else
        disp('WARNING: normalization of fft(pressure) was unsuccessful')
        disp(['t-s: ' sprintf('%.10f',norm(ps))])
        disp(['fft(t-s): ' sprintf('%.10f',norm(pfa))])
        return
        
    end
    psdpf = (1/(fs*n))*abs(pfa).^2;
    psdpf(2:end-1) = 2*psdpf(2:end-1); %psd of ps

    if ~isreal(psurf)
        disp('ifft has returned a complex result')
        psurf = real(psurf);
    end
    if plotfigs
        %plot results
        figure
        plot(psurf,'r'),hold on
        plot(ps)
        xlabel('Samples')
        ylabel('dbar')
        legend('Surface Pressure','Pressure @ Depth')
        prompt = 'Is the adjustment acceptable? [y/n] ';
        result = input(prompt,'s');
        if strcmp(result,'y') || strcmp(result,'yes');
        end
        if strcmp(result,'n') || strcmp(result,'no');
            return
        end
    end
    
    %%4. Calculate Significant Wave Height (H) %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    m0 = psdpf.*freq';
    Hs = 4*sqrt(m0);
    
    %%5. Recalculate h based on attenuation %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    for i = 1:length(psurf)
        h(i,:) = ((((-1.82E-15*psurf(i,:)+2.279E-10)*psurf(i,:)-2.2512E-5)*psurf(i,:)+9.72659)*psurf(i,:))/g;
    end
    %h should not change over the period of averaging
    h = abs(h);h = mean(h);
    
    %%6. Calculate wave group velocity %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %dfn. for Cg from Dalrymple et al. 1984 "Wave Diffraction due to Areas"
    kh = k.*h;
    c = sqrt(g.*tanh(kh))./sqrt(k);
    n = (1/2)*(1+((2*kh)./(sinh(2*kh))));
    Cg = n.*c;
    cgnan = find(isnan(Cg));
    Cg(cgnan) = Cg(cgnan+1);
    
    %%7. Calculate energy (E) %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    E = (1/8)*rho.*g.*Hs.^2;
    
    %%8. Energy flux (F) %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    F = sum(E.*Cg);
    Fp = sum((E.*Cg))./rho;
    Fe = (Fp/sum(Hs.^2));
    
    %%9. Copy data to structure %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    j = 1:length(E);
    WPF.Cg(it,j) = Cg;
    WPF.E(it,j) = E;
    WPF.F(it,:) = F;
    WPF.Fe(it,:) = Fe; %units of epsilon
    WPF.freq(it,j) = freq;
    j = 1:length(psurf);
    WPF.psurf(it,j) = psurf;
    it = it+1;
end
disp('Finished')
WPF.Info.DepStart = [datestr(start,'dd-mm-yyyy HH:MM:SS') ' ICT'];
WPF.Info.DepStop = [datestr(stop,'dd-mm-yyyy HH:MM:SS') ' ICT'];
WPF.Info.avgint = [num2str(intv) ' minute averages'];
WPF.Info.cmt = 'Each row represents a new interval';
WPF.Info.timecmt = 'Times listed are the user-specified start/stop times used to generate these data';
if savedatfile
    %save file
    filename = [name '_' dname];
    sfname = ['Wpf' filename];
    save(sfname,'WPF','-v7.3')
    disp(['Drag force data file saved as ' sfname])
end
% if plotwpf
%     figure
%     [q,~] = size(WPF.F);
%     c = jet(q);
%     for i = 1:q
%         plot(WPF.freq(i,:),WPF.F(i,:),'Color',c(i,:)),hold on
%     end
%     colormap jet
%     cc = colorbar;
%     ylabel('Wave Power Flux (N/s)')
%     xlabel('Frequency (Hz)')
% end
