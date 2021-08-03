function [ub,Tbav]=ubspecfun2(hs,tp,h)
% UBSPECFUN2 - Calculate ub and Tbav from Hs and Tp at surface using JONSWAP
% Same as UBSPECFUN, but with hs,tp,h all matrices (or vectors) of same size
% [ub,Tbav] = ubspecfun2(hs,tp,h)
%
% Input:
%   hs - Significant wave height (m)
%   tp - Peak period (s)
%   h  - Water depth (m)
%
% Written by Pat Wiberg, UVa
% Minor changes by csherwood@usgs.gov
% Vectorization by rsignell@usgs.gov 
g=9.81;
f=30:10:400; f=f./1000; sf=2.*pi.*f; nb=length(f); Tbin=1./f;

[nx,ny]=size(hs);
hs=hs(:);
tp=tp(:);
h=h(:);

% for i=1:length(tp),
%     khd(i)=rtnewt((2*pi./tp(i)).^2*h./g);
% end;
khd=qkhf( 2*pi./tp, h );
khd=khd(:);
ubTd=(pi./tp).*hs./sinh(khd);

sjmat=[];
for i=1:length(hs),
    fp=1./tp(i); sfp=2.*pi*fp;
    gam=1.25; sig=0.08; alp=0.0081;
    if fp~=Inf,
        eterm=-((sf-sfp).^2)/(2.*sig.^2.*sfp.^2);
        ee=exp(eterm);
        t2=gam.^ee;
        t1=-1.25.*(sfp./sf).^4;
        sjmat(i,:)=alp*g.^2./sf.^5.*exp(t1).*t2;
    else,
        sjmat(i,:)=zeros(1,nb);
    end
end;
[np,nbins]=size(sjmat);

kh=zeros(size(sjmat));
for k=1:nb,
    kh(:,k)=qkhf( 2*pi./Tbin(k), h );
    %kh(k)=rtnewt((2*pi./Tbin(k)).^2.*h./g);
end;
ubT=zeros(np,nbins);
HbT=zeros(np,nbins);
waT=zeros(np,nbins);
ub=zeros(np,1);
Hb=zeros(np,1);
Hss=zeros(np,1);
Tbav=zeros(np,1);
fbav=zeros(np,1);

maxiter=10; % max iterations
cf=2*ones(np,1); % initialize cf at 2
for iter=1:maxiter
    ind=find(cf<0.99|cf>1.01); % work only on non-converged points
    if isempty(ind), break, end % quit if nothing left to work on 
    nind=length(ind);
    ubT(ind,:)=(pi./(ones(nind,1)*Tbin).*2.*...
        sqrt(sjmat(ind,:).*0.01)./sinh(kh(ind,:))).^2;
    HbT(ind,:)=4.*sjmat(ind,:).*.01./((cosh(kh(ind,:))).^2);
    waT(ind,:)=4.*sjmat(ind,:).*0.01;
    ub(ind)=sqrt(2).*sqrt(sum(ubT(ind,:),2));  % sum over spectral bands
    Hb(ind)=2.*sqrt(sum(HbT(ind,:),2));
    Hss(ind)=2.*sqrt(sum(waT(ind,:),2));
    
    fbav(ind)=sum(ubT(ind,:)./(ones(nind,1)*(Tbin.^2)),2)./sum(ubT(ind,:),2);
    Tbav(ind)=1./sqrt(fbav(ind));
    cf(ind)=(hs(ind)./Hss(ind));
    sjmat(ind,:)=sjmat(ind,:).*(cf(ind)*ones(1,nbins));
end
ub=reshape(ub,nx,ny);
Tbav=reshape(Tbav,nx,ny);

