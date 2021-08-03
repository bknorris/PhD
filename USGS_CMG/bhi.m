function I = bhi(sv,phi)
% BHI - Calculate big hairy integral in Trowbridge & Elgar, 2001
% I = bhi(sv,phi)
% 
% sv - sigma/V - standard deviation of wave-induced horiz. velocity /
%                magnitude of current velocity 
% phi - angle between waves and currents [radians]
%
% Equation A13 in Trowbridge & Elgar 2001, JP0 31:2402-2417.

% Chris Sherwood, USGS
% February 27, 2002

F=sprintf( '((x.^2-2*(%g)*cos(%g).*x+(%g)).^(1/3).*exp(-0.5 .*x.^2))'...
	   ,1/sv, phi, 1/sv^2);
f = inline(F);
I = quadl(f,-5,5); % not sure this interval is always enough
I = (1./sqrt(2*pi))*(sv).^(2/3).*I;
 
%x = -10:.1:10;
%plot(x,f(x),'-r')