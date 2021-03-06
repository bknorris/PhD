function ncwrite_aqdpWvs(puvData, outFileRoot)
% ncwrite_adcpWvs.m  A function to write Aquadopp wave data, generated by
%                    PUV analysis, to a netCDF file.  
%
%   usage:  ncwrite_aqdpWvs(puvData, outFileRoot);
%       where:  puvData - a structure with the following fields
%                            burst - a vector of burst numbers 
%                               dn - time, ML datenum
%                                p - burst-mean pressure, m
%                             hsig - significant wave height, m
%                            peakF - wave frequency at the peak, Hz
%                          peakDir - wave direction at the peak, deg
%                       peakSpread - wave spreading at the peak, deg
%                               Su - surface elevation spectra based on
%                                    velocity, m^2/Hz
%                               Sp - surface elevation spectra based on
%                                    pressure, m^2/Hz
%                              Dir - wave direction, deg
%                           Spread - wave spreading
%                                F - center frequency of each frequency
%                                    band
%                               dF - bandwidth of each frequency band
%               outFileRoot - a string specifying the name given to the
%                             NetCDF output files, in single quotes
%                             excluding the NetCDF file extension .nc
%
% Written by Charlene Sullivan
% USGS Woods Hole Field Center
% Woods Hole, MA 02543
% csullivan@usgs.gov
%
% Dependencies:
%   julian.m
%   gmean.m
%   gmin.m
%   gmax.m

% C. Sullivan   03/29/06,   version 1.4
% Reverse the chronology on the history attribute so the most recent
% processing step is listed first.  For spectra variables, calculate min
% and max attributes over time for each frequency.
% C. Sullivan   10/28/05,   version 1.3
% ADCP wave direction is direction FROM!!  Wave direction from PUV analysis
% is direction TO, so convert it do direction FROM to be consistent with
% ADCP waves.
% C. Sullivan   06/13/05,   version 1.2
% Wave direction, 'wvdir', is output by PUV analysis as 'direction from'.  
% Convert wave direction to 'direction to' for consistency with ADCP and
% Argonaut peak wave directions.
% C. Sullivan   06/09/05,   version 1.1
% The variables 'pspec' and 'vspec' are converted from m^2/Hz to
% mm/sqrt(Hz)
% C. Sullivan   06/02/05,   version 1.0
% This function writes the data from PUV analysis to the processed data
% NetCDF file.


version = '1.4';

disp(' ')
disp(['Writing statistical wave parameters to ',outFileRoot,'p-cal.nc'])

nc=netcdf([outFileRoot,'p-cal.nc'],'write');

% Convert time from datenum to julian and get
% start_time, stop_time
time_greg = datevec(puvData.dn);
time = julian([time_greg]);
start_time=puvData.dn(1);
stop_time=puvData.dn(end);

% Convert peak wave direction from 'direction to' to 'direction from'
f=find(puvData.peakDir<0);
puvData.peakDir(f)=360+puvData.peakDir(f);
puvData.peakDir(f)=puvData.peakDir(f)-180;

% Convert spectra from m^2/Hz to mm/sqrt(Hz)
puvData.Su = sqrt(puvData.Su') .* 1000;
puvData.Sp = sqrt(puvData.Sp') .* 1000;

% Convert frequency to period
puvData.peakF = 1 ./ puvData.peakF;

% 'frequency' and 'dfreq' are the same for each burst, so just use the
% burst-mean.
puvData.F = gmean(puvData.F');
puvData.dF = gmean(puvData.dF');

% Get the number of records
nRec=length(time);

% Get the number of frequencies
theFreqDim = nc('frequency');
nFreq = size(theFreqDim);
nFreq = nFreq(1);

% Write data to netCDF file
% write coordinate variables
nc{'time'}(1:nRec) = floor(time);                                 
nc{'time2'}(1:nRec) = (time-floor(time)).*(24*3600*1000);   
nc{'burst'}(1:nRec) = puvData.burst;
nc{'lat'}(1)=nc.latitude(:);
nc{'lon'}(1)=nc.longitude(:);

% write record variables
ncvarnames = {'wh_4061','wp_peak','wvdir','spread','hght_18','vspec',...
              'pspec','frequency','dfreq'};
names = {'Hsig','peakF','peakDir','peakSpread','p','Su','Sp','F','dF'};
nNames = length(names);
for i = 1:nNames
    eval(['data = puvData.',names{i},';'])
    eval(['nco = nc{''',ncvarnames{i},'''};'])
    nco.maximum = gmax(data);
    nco.minimum = gmin(data);
    theFillVal = nco.FillValue_(:);
    bads = find(isnan(data));
    data(bads) = theFillVal;
    if strcmp(ncvarnames(i),'vspec') || strcmp(ncvarnames(i),'pspec')
        nco(1:nRec, 1:nFreq) = data; 
    elseif strcmp(ncvarnames(i),'frequency') || strcmp(ncvarnames(i),'dfreq')
        nco(1:nFreq) = data;
    else
        nco(1:nRec) = data;
    end
    disp(['Finished writing ',nco.long_name(:)])
end

% Add the following netcdf file attributes
nc.CREATION_DATE = ncchar(datestr(now,0));
nc.start_time = ncchar(datestr(start_time,0));
nc.stop_time = ncchar(datestr(stop_time,0));

% Update the history attribute
history = nc.history(:);
history_new = ['Data converted to NetCDF by aqdpWvs2nc:ncwrite_aqdpWvs.m V ',...
                version,' on ', datestr(now,0),'; PUV analysis run by ',...
                'aqdpWvs2nc:wad2puv.m V 1.2 on ',datestr(now,0),'; ',history];
nc.history = ncchar(history_new);

% Close netCDF file
nc=close(nc);

disp(['Finished writing statistical wave parameters. ',num2str(toc/60),' minutes elapsed'])

return