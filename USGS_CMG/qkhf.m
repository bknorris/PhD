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

y = (w.^2) .* (h ./G);
% Calculate polynomial on bottom line
kh2 = 1.0 + y.*(D1+y.*(D2+y.*(D3+y.*(D4+y.*(D5+y.*D6)))));

% Calculate final term
kh2 = y.*y + y./kh2;

% return kh
kh = sqrt(kh2);
