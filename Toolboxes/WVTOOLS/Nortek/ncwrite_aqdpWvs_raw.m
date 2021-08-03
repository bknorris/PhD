function ncwrite_aqdpWvs_raw(burstData, outFileRoot)
% ncwrite_adcpWvs_raw.m  A function to write Aquadopp pressure and velocities
%                        to a netCDF file.  
%
%   usage:  ncwrite_aqdpWvs_raw(burstData, outFileRoot);
%       where:  burstData - a structure with the following fields
%                       burst - a vector of burst numbers
%                          dn - time, ML datenum
%                           p - pressure, m
%                           u - east velocity, m/s
%                           v - north velocity, m/s
%                           w - vertical velocity, m/s
%                           a - amplitude, counts
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
%   gmin.m
%   gmax.m

% C. Sullivan   03/29/06,   version 1.2
% Reverse the chronology on the history attribute so the most recent
% processing step is listed first.  For spectra variables, calculate min
% and max attributes over time for each frequency.
% C. Sullivan   10/28/05,   version 1.1
% Provide the user additional feedback regarding code execution.  Add
% version number. Use NetCDF variable objects to write data to the NetCDF
% file.
% C. Sullivan   06/02/05,   version 1.0
% This function writes the raw pressures and velocities to the raw data
% NetCDF file.


version = '1.2';

disp(' ')
disp(['Writing pressure and velocity timeseries to ',outFileRoot,'r-cal.nc'])

nc=netcdf([outFileRoot,'r-cal.nc'],'write');

% Convert time from datenum to julian and get
% start_time, stop_time
nRec=length(burstData.dn(1,:));
for i=1:nRec
    dv = datevec(burstData.dn(:,i));
    time(i,:) = julian([dv]);
end
start_time=burstData.dn(1,1);
stop_time=burstData.dn(end,end);

% Convert velocities from m/s to cm/s for EPIC-compatibility
burstData.u = burstData.u * 100;
burstData.v = burstData.v * 100;
burstData.w = burstData.w * 100;

% Get the number of records
nRec = size(time,1);

% Get the number of samples
theSampleDim = nc('sample');
nSamples = size(theSampleDim);
nSamples = nSamples(1);

% Write data to netCDF file
% write coordinate variables
nc{'time'}(1:nRec, 1:nSamples) = floor(time);                                 
nc{'time2'}(1:nRec, 1:nSamples) = (time-floor(time)).*(24*3600*1000);   
nc{'burst'}(1:nRec) = burstData.burst;
nc{'lat'}(1)=nc.latitude(:);
nc{'lon'}(1)=nc.longitude(:);
nc{'sample'}(1:nSamples) = [1:nSamples]';

% write record variables
ncvarnames = {'hght_18','amp','u_1205','v_1206','w_1204'};
names = {'p','a','u','v','w'};
nNames = length(names);
for i = 1:nNames
    eval(['data = burstData.',names{i},';'])
    eval(['nco = nc{''',ncvarnames{i},'''};'])
    nco.maximum = gmax(data);
    nco.minimum = gmin(data);
    theFillVal = nco.FillValue_(:);
    bads = find(isnan(data));
    data(bads) = theFillVal;
    nco(1:nRec, 1:nSamples) = data;
    disp(['Finished writing ',nco.long_name(:)])
end

% Add the following netcdf file attributes
nc.CREATION_DATE = ncchar(datestr(now,0));
nc.start_time = ncchar(datestr(start_time,0));
nc.stop_time = ncchar(datestr(stop_time,0));

% Update the history attribute
history = nc.history(:);
history_new = ['Pressure and velocity timeseries converted to NetCDF by ',...
               'aqdpWvs2nc:ncwrite_aqdpWvs_raw.m V ',version,' on ',...
               datestr(now,0),'; ',history];
nc.history = ncchar(history_new);

% Close netCDF file
nc=close(nc);

disp(['Finished writing Aquadopp pressure and velocity timeseries. ',...
     num2str(toc/60),' minutes elapsed'])

return
