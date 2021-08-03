%calc wave period
clear
%load any ADV file, in this case, extra pressure and time variables are
%included. Remove extra files and proceed.
load('d:\Projects\Mekong_W2015\Data\Vector\FSS\V5109_070315')
% clearvars -except Pres Time

%%%Basic Settings
start = datenum(2015,03,07,15,30,00);
stop = datenum(2015,03,07,16,00,00);
ind = find(ADV.datetime >= start & ADV.datetime <= stop);

%%%Iterate through time intervals
time = ADV.datetime(ind);pressure = ADV.Pres(ind);
x = (sin(ADV.Metadata.inst_lat/57.29578))^2;
offset = ADV.Metadata.pressure_sensor_height/1000;

window = 120; %120 second averaging interval
step = 60; %60 second step
fs = 32;
avt = step*fs; %samples/step
nwin = window*fs; %samples/window
nsamp = length(time);
ind = [1 avt:avt:nsamp];
T = zeros(length(ind),1);
Ttime = zeros(length(ind),1);
for ii = 1:length(ind)
    if abs(nsamp-ind(ii)) < nwin
        continue
    else
        idx = ind(ii):ind(ii)+nwin-1;
        p = cmgbridge(pressure(idx),100,100,1000);
        ts = time(idx(nwin/2)); %save timestamp of the spectra interval
        %From: UNESCO (1983): Algorithms for computation of fundamental properties
        %of seawater. UNESCO technical papers in marine science 44:1-55.
        g = zeros(length(p),1);h = zeros(length(p),1);
        for i = 1:length(p)
            g(i,:) = 9.780318*(1.0+(5.2788E-3+2.36E-5*x)*x)+1.092E-6*p(i,:);
            h(i,:) = ((((-1.82E-15*p(i,:)+2.279E-10)*p(i,:)-2.2512E-5)*p(i,:)+9.72659)*p(i,:))/g(i,:);
        end
        h = h+offset;
        h = detrend(h);
        
        %Spectra settings
        window = nwin*0.2;
        nfft = window*0.5;
        max_f = 1.2; % maximum frequency
        
        [pf,f] = pwelch(h,hanning(window),round(nfft*0.7),nfft,fs);
        
        ff = find((f<=max_f),1,'last');
        
        pf = pf(2:ff);freq = f(2:ff);
        [~,id] = max(pf);
        if any(pf == 0)
            T(ii,1) = NaN;
        else
            T(ii,1) = 1/freq(id);
        end
        Ttime(ii,1) = time(idx(ii));
    end
end
T(T == 0) = [];Ttime(Ttime == 0) = [];

figure
p(1) = subplot(211);
h(1) = plot(time,pressure,'linewidth',1.5);
set(gca,'xlim',[time(1) time(end)])
legend(h(1),'Pressure')
ylabel('dBar')
datetick('x','dd HH:MM','keepticks','keeplimits')
p(2) = subplot(212);
h(2) = plot(Ttime,T,'color','r','Linewidth',1.5);
set(gca,'xlim',[time(1) time(end)])
legend(h(2),'Peak Period')
ylabel('Seconds')
linkaxes(p,'x')
datetickzoom('x','dd HH:MM','keepticks','keeplimits')