function [Suv,Sp,Dir,Spread,F,dF,DOF] = wavespecdir(u,v,p,fs,zp,zuv,lf,maxfac,minspec,Ndir)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [Suv,Sp,Dir,Spread,F,dF,DOF] = wds(u,v,p,dt,nF,zp,zuv,lf,maxfac,minspec,Ndir);
%
% OUTPUTS:
% Suv is the surface elevation spectra (m^2/Hz), based on the velocity data
% Sp is the surface elevation spectra (m^2/Hz), based on the pressure data
% Dir is the wave direction (deg) 
% Spread is the wave spreading (deg)
% F is the center frequency of each band
% dF is the bandwidth of each band
% DOF is th number of degrees of freedom in each band
%
% INPUTS:
% u,v are detrended east and north components of velocity (m/s), 
% p is pressure (m)
% fs is the sampling frequency (Hz)
% zp is the height of the pressure sensor above the bottom (m)
%     (this means the water depth is the mean pressure plus zp)
% zuv is the height of the velocity cell above the bottom (m)
%
% PARAMETERS:
% lf is the low-frequency cutoff. F<lf are not output
% maxfac is the largest factor scaling pressure to surface elevation
%     spectra and directions at F above this cutoff are returned as NaNs
% minspec is the minimum spectral level for which direction is computed. 
%     directions for spectra below this level are returned as NaNs
% Ndir is the direction of the "north" component (degrees)
% 
% Parameter Recommended values: 
% lf = 0.035;
% maxfac = 200;
% minspec = 0.01;
% Ndir = rotational offset between instrument head and compass; Ndir = 0
% for aquadopps.
%
% Note: The routine is faster if the length of the time series
% factors into small integers, and fastest if it is a power of 2
% 
% Copyright (C) 2001, Lee Gordon, NortekUSA LLC
% Edits made by Ben K Norris, University of Waikato, 2016.
% Changelog: simplified wording (into my own words), hardwired parameters 
% in from original script for clarity about what they do, hardwired nF
% (confusing in original code)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dt = 1/fs;

[np,nt] = size(p);
Dt = np*dt;

f=[1:(np/2)]/Dt;            % frequency array up to Nyquist

uf=(fft(u));             % compute spectrum
vf=(fft(v));             % 
pf=(fft(p));

up=abs(uf(2:(np/2+1),:).^2);       % compute power spectrum & ignore zero freq
vp=abs(vf(2:(np/2+1),:).^2);       % (this uses first half of spectrum)
pp=abs(pf(2:(np/2+1),:).^2);
up=up*2/np^2/f(1);          % scale power spectrum
vp=vp*2/np^2/f(1);
pp=pp*2/np^2/f(1);

pup=real((pf.*conj(uf))*2/np^2/f(1)); % scaled cross-spectra
pup=pup(2:(np/2+1),:);    %limit to the same frequencies as power spectra
pvp=real((pf.*conj(vf))*2/np^2/f(1));
pvp=pvp(2:(np/2+1),:);
puv=real((uf.*conj(vf))*2/np^2/f(1));
puv=puv(2:(np/2+1),:);

nF = length(puv)/10;

[F, Cuu] = logavg(f, up, nF);  % average into log bands
[F, Cvv] = logavg(f, vp, nF);
[F, Cpp] = logavg(f, pp, nF);
[F, Cpu] = logavg(f, pup, nF);
[F, Cpv, dF, Ns, Ne] = logavg(f, pvp, nF);
[F, Cuv] = logavg(f, puv, nF);
DOF=2*(Ne-Ns+1);

aa=find(F>lf);            % low frequency cutoff
lF=length(aa);            % number of frequencies we keep

F=F(aa);
dF=dF(aa);
DOF=DOF(aa);

Cuu=Cuu(aa,:);
Cvv=Cvv(aa,:);
Cpp=Cpp(aa,:);
Cpu=Cpu(aa,:);
Cpv=Cpv(aa,:);
Cuv=Cuv(aa,:);

%find direction from pressure-velocity cross-spectra
Dir = 57.296*atan2(Cpu,Cpv)+Ndir;
id = find(Dir < 0);
Dir(id) = 360+Dir(id);
% Dir=mod(Dir+180,360)-180;

mp=ones(lF,1)*mean(p);        % vertical scaling

F2=F*ones(1,nt);
k=wavek(F2,mp+zp);

sinhkh=sinh(k.*(zp+mp));        %sinh scaling for water depth
coshkh=cosh(k.*(zp+mp));        %cosh scaling for water depth
coshkhz=cosh(k.*zp);            %scaling for pressure sensor elevation
coshkhZ=cosh(k.*(zuv));         %scaling for velocity elevation

ac=find((coshkh.^2)>maxfac);
Suv=(Cuu+Cvv).*(sinhkh./(2*pi*F2)./coshkhZ).^2;
Sp=Cpp.*(coshkh./coshkhz).^2;
ad=union(find(Sp<minspec),ac);

R2=((Cuu - Cvv).^2 + 4*Cuv.^2).^.5./(Cuu+Cvv);
Spread = 57.296*((1-R2)/2).^.5;

% Suv(ac)=nan;
% Sp(ac)=nan;
% Dir(ad)=nan;
% Spread(ad)=nan;



