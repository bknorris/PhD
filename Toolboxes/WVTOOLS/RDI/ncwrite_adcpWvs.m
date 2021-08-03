function ncwrite_adcpWvs(logData,outFileRoot)
% ncwrite_adcpWvs.m  A function to write the timeseries of wave parameters,
%                    output from RD Instruments WavesMon software, to a
%                    netCDF file.  
%
%   usage:  ncwrite_adcpWvs(logData,outFileRoot);
%
%       where:  logData - a structure with the following fields
%                   file    - name of the *_LogData.* file
%                   burst   - burst number
%                   YY      - 2-digit year
%                   MM      - month
%                   DD      - day
%                   hh      - hours
%                   mm      - minutes
%                   ss      - seconds
%                   cc      - 1/100ths seconds
%                   Hs      - significant wave height, meters
%                   Hm      - maximum wave height, meters
%                   Tp      - peak wave period, seconds
%                   Tm      - mean wave period, seconds
%                   Dp      - peak wave direction, degrees true north
%                   ht      - water depth from pressure sensor, millimeters
%               outFileRoot - name given to the processed data NetCDF file
%
% Written by Charlene Sullivan
% USGS Woods Hole Science Center
% Woods Hole, MA 02543
% csullivan@usgs.gov
%
% Dependencies:
%   julian.m
%   gregorian.m
%   gmin.m
%   gmax.m

% C. Sullivan   03/29/06,   version 1.2
% Now using EPIC variable mwh_4064 for maximum wave height. Reverse the
% chronology on the history attribute so the most recent processing step is
% listed first.
% C. Sullivan   10/26/05,   version 1.1
% Provide user additional feedback regarding code execution.  Use NetCDF
% variable objects to handle statistical wave parameters.
% C. Sullivan   06/09/05,   version 1.0
% Do time conversion to julian days in here.  Calculate min/max values and
% replacing nan's with the NetCDF fill value.


version = '1.2';

% Convert time from gregorian to julian and get
% start_time, stop_time, and number of records
time = julian([logData.YY, logData.MM, logData.DD,...
               logData.hh, logData.mm, logData.ss+logData.cc/100]);
start_time=datenum([gregorian(time(1))]);
stop_time=datenum([gregorian(time(end))]);
nrec=length(time);

% Write coordinate variables to NetCDF
nc=netcdf([outFileRoot,'p-cal.nc'],'write');
nc{'time'}(1:nrec) = floor(time);                                 
nc{'time2'}(1:nrec) = (time-floor(time)).*(24*3600*1000);   
nc{'burst'}(1:nrec) = logData.burst;
nc{'lat'}(1) = nc.latitude(:);
nc{'lon'}(1) = nc.longitude(:);

% Write record variables to NetCDF
disp(' ')
disp(['Writing statistical wave parameters to ',outFileRoot,'p-cal.nc'])
ncnames = {'wh_4061','wp_4060','mwh_4064','wp_peak','wvdir','hght_18'};
names = {'Hs','Tm','Hm','Tp','Dp','ht'};
nNames = length(names);
for i = 1:nNames
    eval(['data = logData.',names{i},';'])
    if strcmp(names{i},'ht')
        data = data./1000;
    end
    eval(['nco = nc{''',ncnames{i},'''};'])
    nco.maximum = gmax(data);
    nco.minimum = gmin(data);
    theFillVal = nco.FillValue_(:);
    bads = find(isnan(data));
    data(bads) = theFillVal;
    nco(1:nrec) = data;
    disp(['Finished writing ',nco.long_name(:)])
end

% Add the following NetCDF global attributes
nc.WavesMonCfg.LogDataFile = ncchar(logData.file);
nc.CREATION_DATE = ncchar(datestr(now,0));
nc.start_time = ncchar(datestr(start_time,0));
nc.stop_time = ncchar(datestr(stop_time,0));

% Update the NetCDF history attribute
history = nc.history(:);
history_new = ['Statistical wave parameters converted to NetCDF by ',...
               'adcpWvs2nc:ncwrite_adcpWvs.m V ',version,' on ',...
               datestr(now,0),'; ',history];
nc.history = ncchar(history_new);

% Close NetCDF file
nc = close(nc);

disp(['Finished writing statistical wave parameters. ',...
       num2str(toc/60),' minutes elapsed'])

return
