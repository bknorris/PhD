function e = dissipation_rate(uvw, fs, nfft, iplot) 
% dissipation_rate - Calculate dissipation rate from spectrum
% according to Trowbridge & Elgar 2001, JP0 31:2402-2417.
%
% Input:
%    uvw - 3 x n array of u, v, w velocities
%    fs  - sample rate (Hz)
%    nfft - length of segments used in fft (default = 512)
%
% Output structure e, containing:
%     'e.euv = Dissipation rate calc. from horizontal (u,v) [m^2/s^3]';...
%     'e.eww = Dissipation rate calc. from vertical (w) [m^2/s^3]';...
%     'e.euv_sd = std. dev. of e.euv [m^2/s^3]';...
%     'e.eww_sd = std. dev. of e.eww [m^2/s^3]';...
%     'e.S = mean horizonal current speed [m/s]';...
%     'e.wsd = std. deviation of wave-orbital velocity [m/s]';...
%     'e.uvnoise = noise in U and V [m^2 s^-2 Hz-1]';...
%     'e.IA = Eq. A13 in Trowbridge & Elgar 2001 [ ]';...
%     'e.ss = statistics structure from uvwstats.m'};


% Units in output assume units in input are MKS
if(exist('nfft','var')~=1)
   nfft=512;
   fprintf(1,'Using default nfft=%d/n',nfft);
end
if(exist('iplot','var')~=1)
   iplot=0;
end
alpha = 1.5;    % Kolomogorov constant

ss = uvwstats(uvw(:,1),uvw(:,2));
S = ss.S;
%sd = sqrt(ss.sd1.^2 + ss.sd2.^2);
sd = ss.sd1;
phir = ss.phir;

noverlap = nfft/2;
window = hamming(nfft);
nf = (nfft/2)+1;
[U,f]=pwelch(uvw(:,1),window,noverlap,nfft,fs); U = U(2:nf,1);
[V,f]=pwelch(uvw(:,2),window,noverlap,nfft,fs); V = V(2:nf,1);
[W,f]=pwelch(uvw(:,3),window,noverlap,nfft,fs); W = W(2:nf,1);
f = f(2:nf,1);

% dissipation rate estimates
twopi = 2*pi;
IA = bhi( sd / S, phir );
% Frequency range that is used for inertial-range fit
fm = find(f>=.6 & f<=2.5);
fhn = find(f>=3.5); % assume noise at f > 3.5 Hz
uvnoise  = mean(U(fhn,1)+V(fhn,1));
N = uvnoise*ones(size(U));
Puv = (twopi.*f).^(5/3) .* max(( U+V-N )./twopi, eps*ones(size(U)));
Pww = (twopi.*f).^(5/3) .* W./twopi;
euv = mean( Puv(fm) ) *...
   (55/(21*alpha))*(S)^(-2/3)*(1/IA);
eww = mean( Pww(fm) )*...
   (55/(12*alpha))*(S)^(-2/3)*(1/IA);
euv_sd = std( Puv(fm))*...
   (55/(21*alpha))*(S)^(-2/3)*(1/IA);
eww_sd = std( Pww(fm))*...
   (55/(12*alpha))*(S)^(-2/3)*(1/IA);

e.euv = euv;
e.eww = eww;
e.euv_sd = euv_sd;
e.eww_sd = eww_sd;
e.uvnoise = uvnoise;
e.IA = IA;
e.ss = ss;
e.variables = {...
    'e.euv = Dissipation rate calc. from horizontal (u,v) [m^2/s^3]';...
    'e.eww = Dissipation rate calc. from vertical (w) [m^2/s^3]';...
    'e.euv_sd = std. dev. of e.euv [m^2/s^3]';...
    'e.eww_sd = std. dev. of e.eww [m^2/s^3]';...
    'e.S = mean horizonal current speed [m/s]';...
    'e.wsd = std. deviation of wave-orbital velocity [m/s]';...
    'e.uvnoise = noise in U and V [m^2 s^-2 Hz-1]';...
    'e.IA = Eq. A13 in Trowbridge & Elgar 2001 [ ]';...
    'e.ss = statistics structure from uvwstats.m'};

 
if(iplot)   
   plot(f(fm),euv*ones(size(fm)),'-c','linewidth',3)
   hold on
   plot(f(fm),eww*ones(size(fm)),'-k','linewidth',2)
   plot(f,Puv*(55/(21*alpha))*(S)^(-2/3)*(1/IA),'-b')
   plot(f,Puv*(55/(21*alpha))*(S)^(-2/3),'-b')
   plot(f,Pww*(55/(12*alpha))*(S)^(-2/3)*(1/IA),'color',[.6 .6 .6])
   plot(f,Pww*(55/(12*alpha))*(S)^(-2/3),'color',[.6 .6 .6])
   set(gca,'xscale','log','yscale','log')
end
return