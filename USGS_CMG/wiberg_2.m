function [atype,rtype,E,L]=wiberg_2(Ub_rms,T,d_50)
%
% Function [Type, Eta, Lambda]=ripple1(Ub,T,D)
%
% Function for the prediction of wave-induced
% ripple dimensions for a given wave characteristic
% (Ub=wave orbital velocity; T=wave period)
% and sediment particle size (D in microns)
%
%%%Ub is the significant wave orbital velocity (1.42*sqrt(2)*Ub_rms)

% type='orbital','anorbital'.'suborbital
% 
% If wave height and water depth are known use ripple2.m
%
% Based on the model presented by:
%
% Wiberg, P.L. and Harris, C.K., 1994. Ripple geometry
% in wave-dominated environments, JGR, 99(C1):775-789
%
%
% Ub in m/s
% 
% 

D=d_50*1e-6;  			% Convert particles size from microns to metres
% g=9.81;				% Gravitational acceleration (m/s)
Ub=Ub_rms/100;
do=T.*Ub/pi;      % if orbital velocity is known
%
% Anorbital Ripples L_ano, E_ano
%
L_ano=535*D;
%
% Initial value for do/E_ano=X
%X1=0; 
X0=do./(0.17*L_ano);
X1=zeros(size(X0));
for ij=1:size(X0)
 while abs(X0-X1(ij))>=1e-7
 EoverL(ij,1)=exp(-0.095*(log(X0(ij,1)).^2)+0.442*log(X0(ij,1))-2.28);
 E_ano(ij,1)=EoverL(ij,1)*L_ano;
 X1(ij,1)=do(ij,1)./E_ano(ij,1);
 X0(ij,1)=X1(ij,1);
end
end
%
% Check in do/eta<=10
%
if (do/E_ano)<=10
   E_ano=0.17*L_ano;
end
%
% Criterion for type of ripples
%
L_orb=0.62.*do;
E_orb=0.17.*L_orb;
%
Crit=do/E_ano;
%
if (Crit<20)
%   disp([[' orbital ripples']])
   atype='   orbital';
   rtype=1;
   L=L_orb;
   E=E_orb;
elseif (Crit>100)
%   disp([[' anorbital ripples']])
   atype=' anorbital';
   rtype=3;
   L=L_ano;
   E=E_ano;
elseif (Crit>=20 && Crit<=100)
%   disp([[' suborbital ripples']])
   atype='suborbital';
   rtype=2;
   L=exp( ((log(do/E_ano)-log(100))./(log(20)-log(100)).*(log(L_orb)-log(L_ano))+log(L_ano)) );
   % Initial value for do/E_ano=X
   X1=0; 
   X0=do./(0.17*L);
   while abs(X0-X1)>=1e-7
     EoL=exp(-0.095*(log(X0).^2)+0.442*log(X0)-2.28);
     %E_sub=EoverL*L; %changed by TN
     E_sub=EoL*L;
     X1=do./E_sub;
     X0=X1;
   end
   if (do/E_sub)<=0.17
     E_sub=0.17*L;
  end
  E=E_sub;  
end  