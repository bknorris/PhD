function [u,v] = mag2truenorth(u,v,dec)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Syntax: 
%   
%  [u,v] = mag2truenorth(u,v,dec)
%
%  Inputs: u: beam 1 velocity
%          v: beam 2 velocity
%          dec: magnetic declination to be applied
%
%  Rotate velocitites so they are relative to true N not magnetic North
%
%  Developed by Dr. Julia C Mullarney, University of Waikato, New Zealand 
%  c. 2013
%
%  Additions made by Benjamin K Norris, 2014
%  Edits: cleaned up and organized original code fragments
%         removed extraneous bits of code
%         
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% Rotate to True North coordinates
declination =dec*pi/180; % degrees 2 radians

u2=u;
v2=v;

u= u2.*(ones(size(u2))*cos(declination)) + ...
         v2.*(ones(size(v2))*sin(declination));
v=-u2.*(ones(size(u2))*sin(declination)) + ...
       v2.*(ones(size(v2))*cos(declination));
disp('Magnetic declination applied')
clear u_old v_old
end