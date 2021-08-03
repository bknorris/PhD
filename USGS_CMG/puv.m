function ws =puv(p, u, v, depth, zp, zuv, fs, nfft, rho, window, lf_cutoff, hf_cutoff )
% PUV - Determine wave heights from p, u, v velocity data
% [ws]=puv (p, u, v, depth, zp, zuv, fs, [nfft], [rho], [window], [lf_cutoff], [hf_cutoff] )
%
% Input:
%    p         pressure (dbar) 
%    u         u velocities (m/s)
%    v         v velocities (m/s)
%    depth     Average water depth (m, positive number)
%    zp        Height of pressure sensor off bottom (m, postive number)
%    zuv       Height of veloctity sensor off bottom (m, positve number)
%    fs        Sampling freqency (Hz)
%    nfft      Length of data to window and process [optional; default=512]
%    rho       Water density [optional; default = 1025]
%    window    Data window for spectral calc [optional; default=hanning(nfft)]
%    lf_cutoff Low-frequency cutoff for wave motions (default = 1/50 s)
%    hf_cutoff High-frequency cutoff for wave motions (default = 1/5 s)
%
% Returns:
%    structure ws
%    ws.Hrmsp = Hrms (=Hmo) from pressure
%    ws.Hrmsu = Hrms from u,v
%    ws.ubr = Representative orbital velocity amplitude in freq. band
%             ( lf_cutoff <= f <= hf_cutoff ) (m/s)
%    ws.omegar = Representative orbital velocity (radian frequency)
%    ws.Tr = Representative orbital velocity period (s)
%    ws.phir = Representative orbital velocity direction (angles from
%             x-axis, positive ccw)
%    ws.azr = Represntative orb. velocity direction (deg; geographic
%             azimuth; ambiguous =/- 180 degrees)
%    ws.ublo = ubr in freq. band (f <= lf_cutoff) (m/s)
%    ws.ubhi = ubr in freq. band (f >= hf_cutoff) (m/s)
%    ws.ubig = ubr in infragravity freq. band (lf_cutoff f <= 1/20) (m/s)
%    ws.spread = wave direction spreading (deg)
%    ws.Sp = surface spectra from pressure sensor
%    ws.Suv = surface spectra from velocity

% Chris Sherwood, USGS
% Laura Landerman - correct calc. of phir and azr, June 24, 2003
% Chris Sherwood  - change calc. of wave contributions, May 22, 2004

g = 9.81;
if(exist('rho')~=1),rho = 1025;,end
if(exist('nfft')~=1),nfft = 512;,end
if(exist('window')~=1),window = hanning(nfft);,end
if(exist('lf_cutoff')~=1),lf_cutoff=1/50;,end
if(exist('hf_cutoff')~=1),hf_cutoff=1/5;,end
noverlap=nfft/2;

p = detrend(p);
u = detrend(u);
v = detrend(v);

% compute wave height from velocities
% Determine velocity spectra for u and v
[Gpp, f] = pwelch( rho*g.*p,window,noverlap,nfft,fs);
df = f(3)-f(2);
[Guu, f] = pwelch( u,window,noverlap,nfft,fs);
[Gvv, f] = pwelch( v,window,noverlap,nfft,fs);

% These pairs of time-domain and spectral stats 
% should be very close to equal:
%  varp = var(rho*g*p)
%  varP = sum( Gpp*df )
%  varu = var(u)
%  varU = sum( Guu*df )
%  varv = var(v)
%  varV = sum( Gvv*df )
%  rmsV = sqrt(varV)
%  rmsv = sqrt(mean( v.^2))
%  ubr1 = sqrt(2.*(varu+varv))
%  ubr2 = sqrt(2*(varU+varV))

% determine wave number
omega = 2*pi .* f;
k = qkhf(omega, depth) ./ depth;

%compute linear wave transfer function 
kh = k*depth;
kzp = k*zp;
kzuv = k*zuv;
omega = omega(:);
nf = length(omega);
Hp  = ones(nf,1);
Huv = ones(nf,1);

% change wavenumber at 0 Hz to 1 to avoid divide by zero
i = 1:nf;
if(omega(1)==NaN | omega(1)<=0),
  i = 2:nf;
  Hp(1)=1;
  Huv(1)=1;
end
Hp(i,1)  =      rho*g .* ( cosh(kzp(i,1))  ./ cosh(kh(i,1)) );
Huv(i,1) = omega(i,1) .* ( cosh(kzuv(i,1)) ./ sinh(kh(i,1)) );

% create cut off freqency, so noise is not magnified
ff = min(find(f>=lf_cutoff));
lf = max(find(f<=hf_cutoff));

% combine horizontal velocity spectra
Guv = Guu + Gvv;

% Determine wave height for velocity spectra
Snp = Gpp ./ (Hp.^2);
Snu = Guv ./ (Huv.^2);

% Determine rms wave height (mult by another sqrt(2) for Hs)
% Thornton and Guza say Hrms = sqrt(8 mo)
Hrmsu = 2*sqrt( 2*sum( Snu(ff:lf)*df ) );
Hrmsp = 2*sqrt( 2*sum( Snp(ff:lf)*df ) );

% These are representative orbital velocities for w-c cacluations,
% according to Madsen (1994) Coastal Engineering 1994, Proc., 24th
% Intl. Conf., Coastal Eng. Res. Council / ASCE. pp.384-398.
% (esp. p. 395)
ubr = sqrt( 2*sum(Guv(ff:lf)*df));
omegar = sum( omega(ff:lf).* Guv(ff:lf)*df )./ ...
	 sum( Guv(ff:lf)*df );
Tr = 2*pi/omegar;
% phi is angle wrt to x axis; this assumes Guu is in x direction
% phir = atan2( sum(Guu(ff:lf)*df), sum(Gvv(ff:lf)*df) );  

% this is the line changed on 6/24/03 - I still think it is wrong (CRS)
phir = atan2( sum(Gvv(ff:lf)*df), sum(Guu(ff:lf)*df) );

% convert to degrees; convert to geographic azimuth (0-360, 0=north)
azr = 90-(180/pi).*phir;

% Freq. bands for variance contribs
ig = max( find( f <= 0.05 ));
% low freq, infragravity, high-freq
if(1<ff), 
  ublo = sqrt(2*sum(Guv(1:ff)*df));
else
  ublo = 0;
end
if(ig>ff),
    ubig = sqrt(2*sum(Guv(1:ff)*df));
else
    ubig = 0;
end
if(lf < nfft ),
  ubhi = sqrt(2*sum(Guv(lf:end)*df));
else
  ubhi = 0;
end
% stick results in structure
ws.Hrmsp = Hrmsp;
ws.Hrmsu = Hrmsu;
ws.ubr = ubr;
ws.omegar = omegar;
ws.Tr = Tr;
ws.phir = phir;
ws.azr = azr;
ws.ublo = ublo;
ws.ubhi = ubhi;
ws.ubig = ubig;

if(0),
  figure(1), clf
  subplot(311)
  plot(u)
  subplot(312)
  plot(v)
  subplot(313)
  plot(p)

  figure(2)
  clf
  subplot(121)
  plot(f,Hp)
  ylabel('Hp')
  subplot(122)
  plot(f,Huv)
  ylabel('Huv')
  
  figure(3)
  clf
  plot(f(ff:lf),Snp(ff:lf),'-b')
  hold on
  plot(f(ff:lf),Snu(ff:lf),'-r')
end

function kh = qkhf( w, h )
% QKHF  Quick explicit calculation of kh in dispersion relationship.
%
% kh = qkhf( w, h )
%
% Hard-wired for MKS units.
% Dean and Dalrymple (1991), p 72.
%
% Input:
%  w Angular wave frequency = 2*pi/T where T = wave period [1/s]
%  h Water depth [m]
% Returns:
%  kh = wavenumber * depth [ ]
% 
% Either w or h can be a vector, but not both.

% Chris Sherwood, USGS
% March 17, 1999

D1=0.6666666666;
D2=0.3555555555;
D3=0.1608465608;
D4=0.0632098765;
D5=0.0217540484;
D6=0.0065407983;
G = 9.80665;

y = (w.*w).*h/G;
% Calculate polynomial on bottom line
kh2 = 1.0 + y.*(D1+y.*(D2+y.*(D3+y.*(D4+y.*(D5+y.*D6)))));

% Calculate final term
kh2 = y.*y + y./kh2;

% return kh
kh = sqrt(kh2);







