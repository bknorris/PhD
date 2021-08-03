clear
load('D:\Projects\Mekong_W2015\DataAnalysis\Paper1\F2F3_1Vels.mat')
vp1 = dat.vpro1;clear dat
bin = 5;
win = 180; %seconds (5 minutes)
step = 60; %seconds
avt = 50*step;
nwin = 50*win;
nsamp = length(vp1.time);
ind = [1 avt:avt:nsamp];
for ii = 1:length(ind)
    if abs(nsamp-ind(ii)) < nwin  %skip the last few indexes approaching the end of the t-s
        continue
    else
        idx = ind(ii):ind(ii)+nwin-1;
    end
        %spectra of ADVs
        u = vp1.x(idx,bin);
        v = vp1.y(idx,bin);
        w = (vp1.z1(idx,bin)+vp1.z2(idx,bin))./2;
        time = vp1.time(idx(1));
        
        %calculate RMS velocities (Luhar et al. 2013)
        Ec = (1/nwin)*sum(u);Nc = (1/nwin)*sum(v);
        Ewrms = sqrt((1/nwin)*sum((u-Ec).^2));
        Nwrms = sqrt((1/nwin)*sum((v-Nc).^2));
        Wc = (1/nwin)*sum(w);Wwrms = sqrt((1/nwin)*sum((w-Wc).^2));
        Uc = sqrt(Ec^2+Nc^2);Uwrms = sqrt(Ewrms^2+Nwrms^2);
        
        wave.time(ii) = time;
        wave.Uc(ii) = Uc;
        wave.Uwrms(ii) = Uwrms;
        wave.Wc(ii) = Wc;
        wave.Wwrms(ii) = Wwrms;
end

figure
p(1) = plot(wave.time,wave.Uwrms.*sqrt(2),'r');hold on
p(2) = plot(wave.time,wave.Wwrms.*sqrt(2),'b');
legend(p,{'U_w';'W_w'})
ylabel('Velocity (m/s)')
datetick('x','HH:MM:SS')
xlabel(['Time on ' datestr(wave.time(1),'dd-mm-yy')])
title('Horizontal & vertical wave velocities, h = 0.016 m above bed')