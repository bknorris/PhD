%Plot pressure spectra for depth and surface pressures (non-attenuated and
%attenuated) from HR3, mudflat aquadopp during FSS3.2
%Plot will be used during the ONR Data workshop meeting in HCMC

clear
load HR3_9March2015_f_pad.mat
dname = 'FSS3'; %deployment
name = 'ADHR3'; %instrument
plotfigs = 0;
savefigs = 0;
z = 0.075; %pressure sensor height
intv = 10; %averaging interval (minutes)
start = datenum(2015,03,08,13,45,00);
stop = datenum(2015,03,08,16,25,00);
savefigdir = 'c:\Users\bkn5\Projects\Documents\ONRTalk\';
dat = struct();
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%since the AQDPs recorded over multiple cycles, provide limits for when the
%instrument was underwater.
ind = find(aqdp.datenum >= start & aqdp.datenum <=stop);
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
    
    %Filter timeseries to remove low-f and higher f-noise %%%%%%%%%%%%
    nf = fs/2; %nyquist criterium
    [b,a] = butter(5,0.08/nf,'high');
    p = filtfilt(b,a,p); %filter pressure data
    [b,a] = butter(5,2/nf,'low');
    p = filtfilt(b,a,p); %filter pressure data
    %Compute fft for spectral analysis of pressure data %%%%%%%%%%%%%%
    ps = detrend(p);
    n = length(ps);
    pf = fft(ps,length(p));
    pf = pf(1:n/2+1); %get fft(ps) corresponding to non-negative bandwidths
    freq = 0:fs/length(ps):fs/2;
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
    
    %Calculate wavenumber (k)
    disp('Computing surface pressure attenuation adjustment')
    k = zeros(length(freq),1);
    attn = zeros(length(freq),1); %attenuation
    T = 1./freq; %period = 1/freq
    omega = 2*pi./T';
    for i = 1:length(freq)
        k(i) = om2k(omega(i),h);
        %Calculate pressure attentuation for each f %%%%%%%%%%%%%%%%%%
        %From Raubenheimer, B et al. 1998 "Estimating Wave Heights from
        %Pressure Measured in Sand Bed"
        attn(i) = (cosh(k(i)*(h+z)))/(cosh(k(i)*h));
    end
    attn(attn < 0.2) = 0.2; %eliminate over-reduction less than 0.02Hz
    attn(freq<0.01 | freq>(nf/1.5)) = 1;	% only for selected range (0.01-2.67Hz)
    disp('Applying adjustment to timeseries')
    pfa = pf.*attn; %adjusted fourier transfm.
    psdpfa = (1/(fs*n))*abs(pfa).^2;
    psdpfa(2:end-1) = 2*psdpfa(2:end-1); %psd of psa
    psdpf = (1/(fs*n))*abs(pf).^2;
    psdpf(2:end-1) = 2*psdpf(2:end-1); %psd of ps
    dat.psdpf(it,:) = psdpf;
    dat.psdpfa(it,:) = psdpfa;
    dat.freq = freq;
    it = it+1;
end
t = intv*(it-1);
[num] = max(dat.psdpf(:));
[x,y] = ind2sub(size(dat.psdpf),find(dat.psdpf==num));
integ = trapz(dat.freq,dat.psdpf(x,:));
Hrms = sqrt(integ);
disp(['Hrms = ' num2str(Hrms) 'm'])

%PLOT%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
f1 = figure;
set(f1,'PaperOrientation','portrait',...
    'position',[400 200   1000   800]);
set(gcf,'color','w','PaperPositionMode','auto')
[n,m] = size(dat.psdpf);
c = jet(n);
ax(1) = subplot(2,1,1);
for ii = 1:n
    pq(ii) = plot(dat.freq',dat.psdpf(ii,:),'-x','linewidth',1.5);
    set(pq(ii),'Color',c(ii,:))
    hold on
    box on
%     xlabel('\bf\itHz','FontSize',14)
    ylabel('\bf\itdbar^2/Hz','FontSize',14)
end
grid on
set(gca,'XTickLabel',[])
cb = colorbar('northoutside');colormap jet
caxis([0 t])
set(get(cb,'title'),'string','\bf\itMinutes Elapsed','FontSize',14);
ax(2) = subplot(2,1,2);
for ii = 1:n
    pq(ii) = plot(dat.freq',dat.psdpfa(ii,:),'-x','linewidth',1.5);
    set(pq(ii),'Color',c(ii,:))
    hold on
    box on
    xlabel('\bf\itHz','FontSize',14)
    ylabel('\bf\itdbar^2/Hz','FontSize',14)
end
grid on
set(ax,'Xlim',[0.05 1],'Ylim',[0 0.2],...
    'GridLineStyle',':',...
    'FontName', 'Helvetica','FontSize',14)
set(ax(1),'position',[0.12 0.50 0.8 0.37])
set(ax(2),'position',[0.12 0.1 0.8 0.37])
if savefigs
    fpath = savefigdir;fname = [dname '_' name '_spectra'];
    export_fig([fpath fname],'-png','-native')
    disp(['Figure ' fname '.png saved'])
end

