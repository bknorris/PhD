function [W, Aw, nbw_factor] = hanning_window(N);
%      Generates a hanning window function with N points.  Also returns the
%   PSD correction factor, Aw, and the effective noise bandwidth factor,
%   nbw_factor.  The effective noise bandwidth is df*nbw_factor.
%
%   Usage:  [W, Aw, nbw_factor] = hanning_window(N);
%
B = [1, -1];
args = (2*pi*((1:N) - 1)/N).';
W = B(1) + B(2)*cos(args);
ms = sum(W.*W)/N;
avg = sum(W)/N;
Aw = 1/ms;
nbw_factor = ms/(avg^2);