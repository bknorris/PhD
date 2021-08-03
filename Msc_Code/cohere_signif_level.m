function sl = cohere_signif_level( df, p )
% COHERE_SIGNIF_LEVEL - estimate significance level coherence using
% theoretical formula.
%
% Usage: sig_lev = cohere_bootstrap_signif_level( df, p )
%
% where p is the significance level and df is the number of degrees of
% freedom as described in Shumway and Stoffer (see below). p defaults to 0.95.
%
% sig_lev is the significance level.
%
% This function uses the finv function which is in the Matlab Statistical
% Toolbox.
%
% This function uses an equation taken from Shumway and Stoffer
% (2000), "Time series analysis and its applications", pg. 250,
% equ. 3.82.  
%
% The degrees of freedom is given by:
%
% df = 2Ln/n' when a series of length n has been zero padded to length
% n'.  L is the number of frequencies in the averaging used to compute
% the spectra.  
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% 	$Id: cohere_signif_level.m,v 1.2 2005/01/08 02:49:29 dmk Exp $	
%
% Copyright (C) 2004 David M. Kaplan
% Licence: GPL (Gnu Public License)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin < 2
  p = 0.95;
end

sl = finv( p, 2, df - 2 );
sl = sl / ( df/2 - 1 + sl );
