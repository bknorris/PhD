function theJulian = datenum2julian(theDateNum)

% datenum2julian -- Convert Julian Day to Matlab datenum.
%  datenum2julian(theDateNum) converts theDayNum (Matlab
%   datenum) to its equivalent Julian day.  The Julian
%   day is referenced to midnight, not noon.
 
% Copyright (C) 1998 Dr. Charles R. Denham, ZYDECO.
%  All Rights Reserved.
%   Disclosure without explicit written consent from the
%    copyright owner does not constitute publication.
 
% Version of 26-Oct-1998 15:49:22.

if nargin < 1, help(mfilename), return, end

t0 = datenum(1968, 5, 23) - 2440000;   % May 23, 1968.

result = theDateNum - t0;

if nargout > 0
	theJulian = result;
else
	disp(result)
end
