function [LAM, ETA] = styles_glenn_2002(T,Ub_rms,d_50)
% T in seconds
%Ub in cm/s
%d_50 in microns
nu=0.0119;% Kinematic Viscosity
g=981;
s=2.65;

Ub=sqrt(2)*Ub_rms(1:end,1); %to Ub eq
d_median=d_50.*1e-4; %microns --> cm
Omega=2*pi./T;  % omega=2*pi/T, where T=wave period 
Ab=(Ub)./Omega;

CHI=4*nu*(Ub).^2/(d_median*((s-1)*g*d_median)^(1.5));

Iclt=find(CHI < 2);
Icgt=find(CHI >= 2);

if isempty(Icgt) == 1; %less than 2
  ETA=Ab*0.30.*CHI.^(-0.38);
  LAM=Ab*1.95.*CHI.^(-0.30);
elseif isempty(Iclt) == 1; %greater than 2
  ETA=Ab*0.48.*CHI.^(-1.10);
  LAM=Ab*2.80.*CHI.^(-0.82);
else
  ETA(Iclt)=Ab(Iclt)*0.30.*CHI(Iclt).^(-0.38);
  ETA(Icgt)=Ab(Icgt)*0.48.*CHI(Icgt).^(-1.10);
  LAM(Iclt)=Ab(Iclt)*1.95.*CHI(Iclt).^(-0.30);
  LAM(Icgt)=Ab(Icgt)*2.80.*CHI(Icgt).^(-0.82);
end;

dum=find(isnan(Ub));
ETA(dum)=NaN;
LAM(dum)=NaN;

ETAsg_2002=ETA; %in cm
LAMsg_2002=LAM; %in cm