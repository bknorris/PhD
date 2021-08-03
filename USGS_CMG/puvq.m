function ws =puvq(p, u, v, depth, zp, zuv, fs, nfft, rho, lf_cutoff, hf_cutoff, window )
% PUVQ - Determine wave heights from p, u, v velocity data
% [ws]=puvq (p, u, v, depth, zp, zuv, fs, [nfft], [rho],  [lf_cutoff], [hf_cutoff], [window],)
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
%    ws.Tpp = Peak period from pressure (s)
%    ws.Tpu = Peak period from velocity (s)
%    ws.phir = Representative orbital velocity direction (angles from
%             x-axis, positive ccw)
%    ws.azr = Represntative orb. velocity direction (deg; geographic
%             azimuth; ambiguous =/- 180 degrees)
%    ws.ublo = ubr in freq. band (f <= lf_cutoff) (m/s)
%    ws.ubhi = ubr in freq. band (f >= hf_cutoff) (m/s)
%    ws.ubig = ubr in infragravity freq. band (lf_cutoff f <= 1/20) (m/s) 

% Chris Sherwood, USGS
% Laura Landerman - correct calc. of phir and azr, June 24, 2003
% Chris Sherwood  - change calc. of wave contributions, May 22, 2004
% Chris Sherwood  - calc. params using only ff:lf to avoid divide by zero
%                   msg. June 2, 2006
% Chris Sherwood  - switch to Soulsby method (qkhfs) for wavenumber July 2007
% Chris Sherwood  - add freq. range to computation of omegar July 2007
%                 - embed qkhfs instead of qkhf
%                 - add window arguement to function statement
% Chris Sherwood  - Dec 2008 - revisit direction per J. Lacy suggestions
% Chris Sherwood  - Feb 2008 - add Tpp and Tpu
% Patrick Dickhudt - Nov. 2009 - make Tpp and Tpu nan if any Snp/Snu are
% nan or all are zero - otherwise Tpp/Tpu is max period in window 

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
[Gpp, f] = p_welch( rho*g*p,nfft, fs, window, noverlap );
df = f(3)-f(2);
[Guu, f] = p_welch( u, nfft, fs, window, noverlap);
[Gvv, f] = p_welch( v, nfft, fs, window, noverlap);

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
k = qkhfs(omega, depth) ./ depth;

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
Snp = Gpp(ff:lf) ./ (Hp(ff:lf).^2);
Snu = Guv(ff:lf) ./ (Huv(ff:lf).^2);
fclip = f(ff:lf);

% Determine rms wave height (mult by another sqrt(2) for Hs)
% Thornton and Guza say Hrms = sqrt(8 mo)
Hrmsu = 2*sqrt( 2*sum( Snu*df ) );
Hrmsp = 2*sqrt( 2*sum( Snp*df ) );

% These are representative orbital velocities for w-c cacluations,
% according to Madsen (1994) Coastal Engineering 1994, Proc., 24th
% Intl. Conf., Coastal Eng. Res. Council / ASCE. pp.384-398.
% (esp. p. 395)
ubr = sqrt( 2*sum(Guv(ff:lf)*df));
ubr_check = sqrt(2*var(u)+2*var(v));
omegar = sum( omega(ff:lf).* Guv(ff:lf)*df )./ ...
	 sum( Guv(ff:lf)*df );
Tr = 2*pi/omegar;
% peak frequencies
if any(isnan(Snp)) || ~any(Snp)
    Tpp = nan;
else
    [emax,jpeak]=max(Snp);
    Tpp = 1/fclip(jpeak);
end

if any(isnan(Snu)) || ~any(Snu)
    Tpu = nan;
else
    [emax,jpeak]=max(Snu);
    Tpu = 1/fclip(jpeak);
end

% phi is angle wrt to x axis; this assumes Guu is in x direction
% phir = atan2( sum(Guu(ff:lf)*df), sum(Gvv(ff:lf)*df) );  

% this is the line changed on 6/24/03 - I still think it is wrong (CRS)
% phir = atan2( sum(Gvv(ff:lf)*df), sum(Guu(ff:lf)*df) );

% This is Jessie's replacement for direction
% 12/08 Jessie notes that Madsen uses velocity and suggests
% 12/08 Jessie notes that Madsen uses velocity
%Suu = sqrt(Guu);
%Svv = sqrt(Gvv);
%Suv = sqrt(Guv);
% but I think eqn. 24 is based on u^2, so following is ok:
rr = corrcoef(u,v);
ortest = sign(rr(2,1));
phir = atan2( ortest*sum(Gvv(ff:lf)*df), sum(Guu(ff:lf)*df) );

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
ws.ubr_check = ubr_check;
ws.omegar = omegar;
ws.Tr = Tr;
ws.Tpp = Tpp;
ws.Tpu = Tpu;
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
  plot(f(ff:lf),Hp)
  ylabel('Hp')
  subplot(122)
  plot(f(ff:lf),Huv)
  ylabel('Huv')
  
  figure(3)
  clf
  plot(f(ff:lf),Snp,'-b')
  hold on
  plot(f(ff:lf),Snu,'-r')
end

function kh = qkhfs( w, h )
% QKHFS - Quick iterative calculation of kh in dispersion relationship
% kh = qkhf( w, h )
%
% Input:
%  w Angular wave frequency = 2*pi/T where T = wave period [1/s]
%  h Water depth [m]
% Returns:
%  kh = wavenumber * depth [ ]
% 
% Either w or h can be a vector, but not both.
% Hard-wired for MKS units.
% Orbital velocities from kh are accurate to 3e-12 !
%
% RL Soulsby (2006) "Simplified calculation of wave orbital velocities"
% HR Wallingford Report TR 155, February 2006
% Eqns. 12a - 14

% csherwood@usgs.gov
% Sept 10, 2006

g = 9.81;
x = w.^2*h./g;
y = sqrt(x) .* (x<1) + x.* (x>=1);
%this appalling bit of code is about 25% faster than a for loop
t = tanh( y );
y = y-( (y.*t -x)./(t+y.*(1-t.^2)));
t = tanh( y );
y = y-( (y.*t -x)./(t+y.*(1-t.^2)));
t = tanh( y );
y = y-( (y.*t -x)./(t+y.*(1-t.^2)));
kh=y;
return;

