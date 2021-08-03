function wstats = wave_stats(spect,freq,max_f,varargin)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Script compiled by B K Norris c. 2015, adapted from Urs Neumeier, 2003
%version 1.0
%
%
% Usage: wstats = wave_stats(spect,freq,max_f,varargin)
%
% INPUT:
% spect: power spectral density of pressure signal
% freq:  frequency range given for spectrum (0:nf)
% max_f: upper frequency limit
% 
% If you want zero_crossing parameters, enter the following values:
% 
% p = sea surface elevation time series (detrended)
% h = mean water depth (m)
% fs = sampling frequency (Hz)
% z = height of pressure sensor from seabed (m)
% noseg = number of segments for fft (optional, default: 4)
% Corr_lim = [min max] frequency for attenuation correction (Hz, 
%            optional, default [0.05 0.33])
%
% OUTPUT:
% wstats: contains wave statistics (Hm0, Hs, etc.)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% narginchk(3,9)
if length(varargin) == 0
    do_zero_crossing = 0;
end
if length(varargin) ~= 0 
    if length(varargin) < 3
        error('Not enough input values for zero crossing parameters')
    else
        p = varargin{1};
        h = varargin{2};
        fs = varargin{3};
        do_zero_crossing = 1; %turn on flag for zero-x calcs
    end
    if length(varargin) < 4
        warning('No value for z is provided')
        prompt = 'Enter a value for z (in m): ';
        z = input(prompt);
    else
        z = varargin{4};
    end
    if length(varargin) < 5 || isempty (varargin{5})
        noseg = 4;	    			% No of segments (for spectrum)
    else
        noseg = varargin{5};
    end
    if length(varargin)< 6 || isempty (varargin{6})% process Corr_lim or take default values
        min_f = 0.05;		% mininum frequency, below which no correction is applied (0.05)
        max_f = 0.33;		% maximum frequency, above which no correction is applied (0.33)
    else
        Corr_lim=varargin{6};
        if diff(Corr_lim)<=0 | any(Corr_lim<0) | length(Corr_lim)~=2
            error('Incorrect Corr_lim argument');
        end
        min_f = Corr_lim(1);
        max_f = Corr_lim(2);
    end
end
    
if ~exist('h','var')
    h = 0;
end
% frequence range over which the spectrum is integrated for calculation of the moments
integmin=min(find(freq >= 0));	    % this influences Hm0 and other wave parameters
integmax=max(find(freq <= max_f*1.5));

df = freq(1);						% bandwidth (Hz)
for i=-2:4							% calculation of moments of spectrum
	moment(i+3)=sum(freq(integmin:integmax).^i.*spect(integmin:integmax))*df;
end
         
m0 = moment(0+3);
Hm0 = 4 * sqrt(m0);                 % Hsig by spectral means
[~,E] = max(spect);                 % Frequency at max of spectrum
fpeak = freq(E);											
Tp = 1/fpeak;                       % Peak period

%*************************************************************
% Calculation of zero-crossing parameters
            % length of segment
if exist('pr_corr.m','file') && exist('zero_crossing.m','file') && do_zero_crossing
    m = length(p);				
    n = fix(m/noseg/2)*2;	
	if ~isempty(z)				    % if z is empty, no attenuation correction
 	   pt_surf = pr_corr(p,h,fs,z,n,[min_f max_f]); %correct p for attenuation
	else
 	   pt_surf=p;
	end
	[resZcross,namesZcross]=zero_crossing(pt_surf,fs); % do the zero-crossing analysis
else
	resZcross=[];
	namesZcross={};
end

%*************************************************************
% Output the results

% put results together
res = [h,Hm0,Tp,m0,resZcross];
names = {'h','Hm0','Tp','m0',namesZcross{:}};
for i = 1:length(names)
    wstats.(names{i}) = res(i);
end
% disp('Displaying Wave Statistics')
% for i=1:length(res)
%     fprintf('%14s:   %g\n',names{i},res(i));
% end

end