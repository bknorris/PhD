function [v1,varargout] = adjsound(temp,sal,ss,v1,varargin)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%	Adjust acoustic instrument velocities based on a temperature salinity
%	relationship for the speed of sound.
%
%   temp: a vector of temperature measurements
%   sal: a vector of salinity measurements
%   ss: a vector of instrument soundspeed
%   v1-v4: velocity measurements (beams 1-4), at least one velocity
%   measurement is required
%
%   Usage: [v1,v2,v3,v4] = adjsound(temp,sal,ss,v1,v2,v3,v4)
%
%   This script was written by Benjamin K Norris, 2015
%   University of Waikato, New Zealand
%  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Set up the defaults
if nargin<4
    help(mfilename)
    return
end

t = temp./10;
S = sal;
c = 1449.05 + 45.7.*t - 5.21.*t.^2 + 0.23.*t.^3 + (1.333 - 0.126.*t + 0.009.*t.^2).*(S - 35);
%t = T/10 where T = temperature in degrees Celsius, D = 0;
%From:
%A.B. Coppens, Simple equations for the speed of sound in Neptunian waters (1981)
%J. Acoust. Soc. Am. 69(3), pp 862-863
ssnew = mean(c);
ssold = mean(ss);
v1 = v1.*(ssnew/ssold);
if length(varargin) == 1; %#ok<*ISMT>
    v2 = varargin{1};
    v = v2.*(ssnew/ssold);
    for k = 1:length(varargin)
        varargout{k} = v(:,k); %#ok<*AGROW>
    end
elseif length(varargin) == 2;
    v2 = varargin{1};
    v3 = varargin{2};
    v2 = v2.*(ssnew/ssold);
    v3 = v3.*(ssnew/ssold);
    varargout{1} = v2;
    varargout{2} = v3;
elseif length(varargin) == 3;
    v2 = varargin{1};
    v3 = varargin{2};
    v4 = varargin{3};
    v2 = v2.*(ssnew/ssold);
    v3 = v3.*(ssnew/ssold);
    v4 = v4.*(ssnew/ssold);
    varargout{1} = v2;
    varargout{2} = v3;
    varargout{3} = v4;
end
disp('T-S correction for velocities applied')
end