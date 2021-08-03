function xx = fill_2D_array(x, Nfft, W, Noverlap)
%
%   Construct a 2-D array from the input time series, x, 
% in anticipation of spectrum processing.  The 2-D array will
% have time records of Nfft points per record with one record
% per column.  A window function can be applied to each record
% if desired.  The input vector, x, is not zero-filled to
% make the last record full-length; a partially filled last record
% is dropped.
%
%   Noverlap is the (optional) overlap factor:  50% overlap 
% is Noverlap = 2; 75% overlap is Noverlap = 4; 90% overlap is 
% Noverlap = 10, etc.  W is the (optional) window function. 
%
%  Usage:   xx = fill_2D_array(x, Nfft);
%     for no overlap and no window
%
%  Usage:   xx = fill_2D_array(x, Nfft, W);
%     for window but no overlap
%
%  Usage:   xx = fill_2D_array(x, Nfft, W, Noverlap);
%     for overlap and window
%
[Nrow, Ncol] = size(x);
if Ncol ~= 1;
    x = x.'; Nt = Ncol;
else
    Nt = Nrow;
end
%
%   Use no overlap if less than 4 input arguments
if nargin < 4
    Noverlap = 1;
end
%
%   Use no window if less than 3 input arguments
if nargin < 3
    W = 1;
end
%
[Nrow, Ncol] = size(W);
if Ncol ~= 1;
    W = W.';
end
%
%  Calculate the shift in number of points from record to record.
% If the result of Nfft/Noverlap is not an integer, the overlap
% will not be exactly the fraction requested.  Nshift must be an
% integer so the 'round' operation is necessary.
Nshift = round(Nfft/Noverlap);
%  Find the total number of records.  (We'll just throw away
% the last record if it is incomplete.)
Nrecs = floor(Nt/Nshift) - (Noverlap-1);
%
%  Initialize the indexes
n1 = 1; n2 = Nfft;
%  Allocate space and define the size of the 2D time-record array
xx = zeros(Nfft,Nrecs);
%
%  Fill time-record array (and apply time-domain window, W)
for ii = 1:Nrecs
    xx(:,ii) = x(n1:n2).*W;
    n1 = n1+Nshift; n2 = n2+Nshift;
end
%
