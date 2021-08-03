%Surface Tracking for Aquadopps

function aqdp = aqdpstrack(a,lat)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Syntax: 
%   
%      aqdp = aqdpstrack(aqdp)
%
%      Surface Tracking for Aquadopp data
%      INPUTS: aqdp - a structure containing processed aquadopp data
%              lat - latitude of instrument (in dd)
%
%      This script uses the pressure signal from the aquadopp to clip all
%      above water bins. It only trims the velocity data, not correlations or
%      backscatter. Depth calculations are based on 'Algorithms for 
%      computation of fundamental properties of seawater, 1983 Unesco'
%
%  Script developed by Benjamin K Norris, University of Waikato, NZ 2014
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
prs = a.pressure;
x = (sin(lat/57.29578))^2;
g = 9.780318*(1.0+(5.2788E-3+2.36E-5*x)*x)+1.092E-6.*prs;
depth = ((((-1.82E-15.*prs+2.279E-10).*prs-2.2512E-5).*prs+9.72659).*prs)./g;
%depth in m

dhts = repmat(a.rangebins,length(depth),1);
bursts = 1:length(a.u);

for ii = 1:length(depth)
    indx = dhts(ii,:) >= depth(ii);
    a.u(ii,indx) = NaN;
    a.v(ii,indx) = NaN;
    a.w(ii,indx) = NaN;
    count = find(indx == 1);
end
aqdp = a;
disp('Surface Tracking Complete')
end