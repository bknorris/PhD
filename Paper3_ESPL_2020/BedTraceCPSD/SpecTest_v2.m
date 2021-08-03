%try filtering the bd level first, then computing spectra.

clear
%load just one of the VPs from the experiment
datdir = 'f:\Mekong_W2015\DataAnalysis\Paper3\VPs\09-03-15\';
if ~exist([datdir 'March9_VP2_XY_bd.mat'],'file')
    load('D:\Projects\Mekong_W2015\Data\Vectrino\9March2015\9March2015_VP2.mat');
    bds = VPRO.Data.BottomCheck_BottomDistance;
    bdt = VPRO.Data.BottomCheck_HostTimeMatlab;
    vpt = VPRO.Data.Profiles_HostTimeMatlab;
    [VPRO.Data,VPRO.Config] = VPQC(VPRO.Data,VPRO.Config,1,1,0); %warning this step is SLOW
    [VPRO.Data,VPRO.Config] = VPro_coordinateTransform(VPRO.Data,VPRO.Config,'bx');
    vpx = VPRO.Data.Profiles_VelX;
    vpy = VPRO.Data.Profiles_VelY;
    data.bdt = bdt;
    data.bds = bds;
    data.vpt = vpt;
    data.x = vpx;
    data.y = vpy;
    clear VPRO
    save([datdir 'March9_VP2_XY_bd'],'-struct','data')
else
    load([datdir 'March9_VP2_XY_bd.mat'])
end
%% Begin Program
bins = 13:17;
start = datenum(2015,03,09,2,10,00);
stop = datenum(2015,03,09,7,40,00);
%rotate velocities to x-shore direction
heading = -20;th = heading*pi/180;
R = [cos(th) -sin(th); sin(th) cos(th)];
vpx = nanmean(x(:,bins),2);
vpy = nanmean(y(:,bins),2);
vpxy = [vpx vpy];rxy = zeros(size(vpxy));
for jj = 1:length(vpxy)
    rxy(jj,:) = vpxy(jj,:)*R;end %rotate x, y to cross & along-shore

%remove NaNs, gap fill with zeros
tid = find(vpt>=start&vpt<=stop);
rxy = rxy(tid,:);
vpt = vpt(tid,:);
tid = find(bdt>=start&bdt<=stop);
bds = bds(tid);
bdt = bdt(tid);
%Downsample velocity to 10Hz
[vpt,idx]=unique(vpt);
x10 = interp1(vpt,rxy(idx,1),bdt);nid = find(isnan(x10));
disp(['Found ' num2str(length(nid)) ' NaNs in x10'])
x10(nid) = 0;
y10 = interp1(vpt,rxy(idx,2),bdt);nid = find(isnan(y10));
disp(['Found ' num2str(length(nid)) ' NaNs in y10'])
y10(nid) = 0;
bds = cmgbridge(bds,1E5,1E5,1E6);
nid = find(isnan(bds));
disp(['Found ' num2str(length(nid)) ' NaNs in bottom trace'])
bds(nid) = 0;
%% Spectral analysis
if mod(length(x10),2)==1 %make ts even
    x10 = x10(1:end-1);
    y10 = y10(1:end-1);
    bds = bds(1:end-1);
end
[np,nt]=size(x10);
x10 = detrend(x10);
y10 = detrend(y10);
bds = detrend(bds);

nF = 200;
fs = 10;
dt = 1/fs;
Dt=np*dt;
f=[1:(np/2)]/Dt;         % frequency array up to Nyquist
w = 0.54-0.46*cos(2*pi*[0:np-1]'/(np-1)); %hamming window

uf=(fft(x10.*w));             % compute spectrum
vf=(fft(y10.*w));
bf=(fft(bds.*w));

up=abs(uf(2:(np/2+1),:).^2);       % compute power spectrum & ignore zero freq
vp=abs(vf(2:(np/2+1),:).^2);       % (this uses first half of spectrum)
bp=abs(bf(2:(np/2+1),:).^2);
up=up*2/np^2/f(1);          % scale power spectrum
vp=vp*2/np^2/f(1);
bp=bp*2/np^2/f(1);

pub=real((bf.*conj(uf))*2/np^2/f(1)); % scaled cross-spectra of x and bd
pub=pub(2:(np/2+1),:);    %limit to the same frequencies as power spectra

msc = ((pub).^2)./(up.*bp);

[~, Cuu] = logavg(f, up, nF);  % average into log bands
[~, Cvv] = logavg(f, vp, nF);
[~, Cbb] = logavg(f, bp, nF);
[~, Msc] = logavg(f, msc, nF);
[F, Cub, Ns, Ne] = logavg(f, pub, nF);
id = find(F<1E-3);
F(id) = NaN;
Cuu(id) = NaN;
Cbb(id) = NaN;


figure(1)

b(1) = loglog(F,Cuu,'k');hold on
b(2) = loglog(F,Cbb,'r');
set(b,'linewidth',1.5)
ylabel('Spectral density')
xlabel('f (Hz)')
legend(b,{'X-shore velocity';'BLE'})
set(gca,'xlim',[10^-4 10^1])
prettyfigures('text',13,'labels',14,'box',1)

DOF=2*(Ne-Ns+1);
v = DOF;
alfa = 1-0.95;
ci = zeros(length(DOF),2);
for j = 1:length(DOF)
    c = chi2inv([1-alfa/2 alfa/2],v(j));
    ci(j,:) = v(j)./c;
end
s2est = std(x10);
ff = sum(w.*w)/np;
cs = (ff*s2est)/(2*(np/2+1)*(1/fs));
ci = ci.*cs;
cf = [Cuu-ci(:,1) Cuu+ci(:,2)];cf(cf<0) = NaN;
nid = ~isnan(cf(:,1));
L = mean(10*log10(Cuu(nid))-10*log10(cf(nid,1)));
U = mean(10*log10(cf(nid,2))-10*log10(Cuu(nid)));
er = ploterr(2,5*10^-2,(5*10^-2)*log10(L),(5*10^-2)*log10(U),'.','logxy','hhxy',0.1);
set(er,'linewidth',1.5,'color','k')
export_fig(['g:\GradSchool\DataAnalysis\Paper3\WorkingFigures\Results\' 'U_bd_spectra_030915'],'-pdf')
%



%run spectra
% M = 600*10;
% N = length(bds);
% sw = hamming(M);
% fs = 10;
% [s12,f] = cpsd(x10,bds,sw,M*0.5,M,fs);
% [s1,~] = pwelch(x10,sw,M*0.5,M,fs);
% [s2,~] = pwelch(bds,sw,M*0.5,M,fs);
% [msc,~] = mscohere(x10,bds,sw,M*0.5,M,fs);
% coh = abs(s12).^2 ./ (s1.*s2);
% phase= atan2(-imag(s12),real(s12))*180/pi;
% nu = 2*floor(N/M); %DOF
% ci = cohere_signif_level(nu);



% figure(1)
% sp(1) = subplot(211);
% loglog(f,s1,'-k')
% title('Velocity')
% ylabel('Spectral density [m^2s^{-2}Hz^{-1}]')
% xlabel('f (Hz)')
% sp(2) = subplot(212);
% loglog(f,s2,'-k')
% set(sp,'xlim',[10^-2 10^1])
% title('Bed Level Elevation')
% ylabel('Spectral density [m^2Hz^{-1}]')
% xlabel('f (Hz)')
%
% figure(2)
% sp(1) = subplot(211);
% semilogx(f,msc,'-k'), hold on
% plot(linspace(10^-2,10^1,10),ones(10,1)*ci,'-r')
% ylabel('Coherence^2')
% xlabel('f (Hz)')
% sp(2) = subplot(212);
% semilogx(f,phase,'-k'), hold on
% plot(linspace(10^-2,10^1,10),zeros(10,1),'-r')
% ylabel('Phase (deg)')
% xlabel('f (Hz)')
% set(sp,'xlim',[10^-2 10^1])
