%Test bed for spectral analysis, bottom height and orbital velocity
clear
close all
load('/Volumes/MEKONG3/Mekong_W2015/DataAnalysis/Paper3/BottomTrack/VTA1a_bdtrace.mat')
bd = vpro1.bdist;bdt = vpro1.time;clear vpro1 vpro2 vpro3
load('/Volumes/MEKONG3/Mekong_W2015/DataAnalysis/Paper3/WaveOrbitalVs/VTA1a_uorb.mat')
ub = vpro1.ubr;vt = vpro1.time;clear vpro1

t1 = vt(1);t2 = vt(end);
id = find(bdt>= t1 & bdt <= t2);
ubr = ub(id);vt = vt(id);
bd = bd(id);bdh = 0.07-bd;t = bdt(id);
bds = fastsmooth(bdh,floor(length(bdh)/20),1,1);
brsd = bsxfun(@minus,(bdh-bds),mean(bdh-bds)); %mean removal & detrend
ubs = fastsmooth(ubr,floor(length(ubr)/20),1,1);
ursd = detrend(ubr-ubs); %mean removal & detrend

figure
p(1) = subplot(411);
plot(t,bdh,'k','linewidth',1.5),hold on
plot(t,bds,'r')
p(2) = subplot(412);
plot(t,brsd,'b')
p(3) = subplot(413);
plot(t,ubr,'k','linewidth',1.5),hold on
plot(t,ubs,'r')
p(4) = subplot(414);
plot(t,ursd,'b')
set([p(1) p(2) p(3)],'xticklabel','')
datetickzoom('x','HH:MM','keepticks','keeplimits')
linkaxes(p,'x')

%Use a sample of the data
b = brsd(20*60*10:40*60*10);
u = ursd(20*60*10:40*60*10);
ts = t(20*60*10:40*60*10);
figure
clear p
p(1) = subplot(211);
plot(ts,b,'k')
title('Bottom Trace')
ylabel('Bed elevation (m)')
p(2) = subplot(212);
plot(ts,u,'k')
title('Wave orbital velocity')
ylabel('m/s')
xlabel(['Time on ' datestr(ts(1)),'dd-mm-yy'])
set(p(1),'xticklabel',[])
datetickzoom('x','HH:MM','keepticks','keeplimits')
linkaxes(p,'x')

%Check for stationarity!
b = diff(b,2);u = diff(u,4);
figure
clear p
p(1) = subplot(211);
[acor,lag]=xcorr(b,'coeff');
l = floor(length(acor)/2)+1;
stem(lag(l:end),acor(l:end)),hold on
vcrit = sqrt(2)*erfinv(0.95); %95 ci
lconf = -vcrit/sqrt(length(b));
upconf = vcrit/sqrt(length(b));
plot([0 length(acor)],[1 1]'*[lconf upconf],'r')
title('Autocorrelation after differencing, bottom trace')
xlabel('k (Lags)')
ylabel('r_k')
p(2) = subplot(212);
[acor,lag]=xcorr(u,'coeff');
l = floor(length(acor)/2);
stem(lag(l:end),acor(l:end)),hold on
plot([0 length(acor)],[1 1]'*[lconf upconf],'r')
set(p,'xlim',[-0.5 20])
linkaxes(p,'xy')
title('Autocorrelation after differencing, u_b_r')
xlabel('k (Lags)')
ylabel('r_k')

%try computing the spectra
b = b(1:length(u));
M = hanning(1000);
dof = 2*floor(length(b)/1000);
[pb,f]=pwelch(b,M,[],6000,10);
[pe,~]=pwelch(u,M,[],6000,10);
[cbe,~]=cpsd(b,u,M,[],6000,10);
[msc,~]=mscohere(b,u,M,[],6000,10);

figure
title('PSD of bottom trace and u_b_r')
xlabel('f (Hz)')
ylabel('PSD')
s(1) = loglog(f,pb);hold on
s(2) = loglog(f,pe,'r');
alfa = 1 - 0.95;
c = chi2inv([1-alfa/2 alfa/2],dof);
c = dof./c; %spectra confidence intervals
cx = c(2)-c(1);
sc = 5^-15;
plot(0.1,10^-10,'.k')
plot([0.1 0.1],[10^-10-(cx/2)*sc 10^-10+(cx/2)*sc],'-k')
text(0.12,10^-10,'95%')
lg = legend(s,{'bottom trace';'u_b_r'});
set(lg,'position',[0.7 0.2 0.05 0.05])
set(gca,'xlim',[10^-3 1.2])
figure
semilogx(f,msc),hold on
ci = cohere_signif_level(30);
plot(f,ones(length(f),1)*ci,':',...
    'Color','k')
xlabel('f (Hz)')
ylabel('Coherence')
title('MSC of bottom trace & u_b_r')
close(figure(1)),prettyfigures('text',12,'labels',12,'box',1)
tilefigs
sdir = '/Volumes/MEKONG3/Mekong_W2015/Figures/Paper3/BedCPSD/';
fname = 'VTA1a';
h = findobj('type','figure');
% export_fig(h(1),[sdir fname 'MSC_bd_ubr'],'-png')
% export_fig(h(2),[sdir fname 'PSD_bd_ubr'],'-png')
% export_fig(h(3),[sdir fname 'acor_bd_ubr'],'-png')
% export_fig(h(4),[sdir fname 'bd_ubr_ts'],'-png')