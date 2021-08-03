function psd = Gxx_2D(array, fs, Aw)
%     Computes the single-sided PSD, psd, from a 2D matrix of time series,
%   array, sampled at fs samples per second.  If the time-series data in
%   'array' were windowed, you must supply the third argument, Aw, which is
%   the psd correction factor for the window function used.
%
%   Usage:  psd = Gxx_2D(array, fs);     % If no time-domain window used
%           psd = Gxx_2D(array, fs, Aw); % If a time-domain window was used
%
if nargin < 3; Aw = 1; end
%
[N, Nrecs] = size(array);
even = true;
if rem(N, 2) ~= 0; even = false; end
%
dt = 1/fs;  T = N*dt;  df = 1/T;
%
lsp = fft(array)*dt;
%
if even
    Nhalf = N/2 + 1;
else
    Nhalf = (N + 1)/2;
end
%
lsp_half = lsp(1:Nhalf,:);
clear lsp
%
psd = 2*abs(lsp_half).^2/T;
psd(1, :) = psd(1,:)/2;
%
if even
    psd(Nhalf, :) = psd(Nhalf, :)/2;
end
%
psd = psd*Aw;