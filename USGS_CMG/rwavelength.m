function rl = rwavelength( Tp, Hd, Vr, phi );
% RWAVELENGTH - Relative wavelength of shallow-water wave with current
% rl = rwavelength( Tp, Hd, Vr, phi );
%
% Tp = period (s)
% Hd = water depth (m)
% Vr = depth-mean current speed (m/s)
% phi = angle from current direction to wave origin (CCW; radians)
%
% Calls Matlab FZERO routine

% Chris Sherwood, USGS
% Last revised 11/17/2003

global Tp_r Vr_r phi_r Hd_r
Tp_r = Tp; Vr_r = Vr; phi_r = phi; Hd_r = Hd; 
Lguess = (2*pi)./ ( qkhf( 2*pi/Tp_r, Hd_r )./Hd_r );
rl = fzero(@rwaveno_func,Lguess);
clear GLOBAL Tp_r Vr_r phi_r Hd_r

function y = rwaveno_func( L )
% RWAVENO_FUNC - Relative wave number function for use with FZERO
% y = rwaveno( L )
% L ( m ) is guesstimated wavelength. y = 0 when correct L is guessed.
%
% Example of calling function:
%
% global Tp_r Vr_r phi_r Hd_r
% Tp_r = 10; Vr_r = 0.2; phi_r = 60*(pi/180); Hd_r = 9; 
% Lguess = (2*pi)./ ( qkhf( 2*pi/Tp_r, Hd_r )./Hd_r );
% L = fzero(@rwaveno,Lguess)
%
% where Tp_r = period (s), Vr_r = depth-mean velocity (m/s),
%       phi_r = angle CCW from current to waves (radians), Hd_r = depth
%       (m)
% Chris Sherwood, USGS
global Tp_r Vr_r phi_r Hd_r
g = 9.80665;
twopi = 2*pi;
y= (L ./Tp_r - Vr_r*cos(phi_r))^2 - (g*L./twopi)*tanh(twopi*Hd_r/L);