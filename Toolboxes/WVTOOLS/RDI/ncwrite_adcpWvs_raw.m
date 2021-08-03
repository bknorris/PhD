function ncwrite_adcpWvs_raw(burstData, outFileRoot);
% ncwrite_adcpWvs_spec.m  A function to write the timeseries of raw data,
%                         output from RD Instruments WavesMon v. 2.01, to
%                         a netCDF file.  
%
%   usage:  ncwrite_adcpWvs_raw(burstData, outFileRoot);
%
%       where:  burstData - an structure w/ the following fields:
%                   burstData.press.data - the raw pressure timeseries, mm
%                   burstData.press.time - time associated with pressure,
%                                          YYYYMMDDhhmmsscc
%                   burstData.strk.data - the raw range-to-surface track
%                                         timeseries, mm
%                   burstData.strk.time - time associated with range-to-
%                                         surface track, YYYYMMDDhhmmsscc
%                   burstData.vel.data - the raw velocity timeseries (mm/s)
%                   burstData.vel.time - time associated with velocity, 
%                                        YYYYMMDDhhmmsscc
%               outFileRoot - name given to the raw data NetCDF file
%
% Written by Charlene Sullivan
% USGS Woods Hole Science Center
% Woods Hole, MA 02543
% csullivan@usgs.gov
%
% Dependencies:
%   julian.m
%   gmin.m
%   gmax.m

% C. Sullivan   03/29/06,   version 1.2
% Reverse the chronology on the history attribute so the most recent
% processing step is listed first. In calculating max and min attributes,
% calculate the max and min over time for each frequency bin.
% C. Sullivan   10/26/05,   version 1.1
% Provide the user more feedback regarding code execution.
% C. Sullivan   06/09/05,   version 1.0
% Do time conversion to julian days in here.  Calculate min/max values. Do
% nothing with respect to replacing bad data values with fill values.
% This means the data in this raw data netcdf file is an exact
% representation of the data in the raw data ascii output files.


version = '1.2';

disp(' ')
disp(['Writing pressure and velocity timeseries to',...
       outFileRoot,'r-cal.nc'])

nc=netcdf([outFileRoot,'r-cal.nc'],'write');

% Get dimensions from NetCDF file
theRecDim = nc('time');
nRec = size(theRecDim);
nRec = nRec(1);
theSampleDim = nc('sample');
nSamples = size(theSampleDim);
nSamples = nSamples(1);
theBeamDim = nc('beam');
nBeams = size(theBeamDim);
nBeams = nBeams(1);
theBeamBinDim = nc('beambin');
nBeamBins = size(theBeamBinDim);
nBeamBins = nBeamBins(1);
clear *Dim

% Check raw data times, which ought to be equal
time_check = isequal(burstData.press.time, burstData.strk.time, burstData.vel.time);
if true(time_check)
    raw_time = burstData.press.time';
    time = YYYYMMDDhhmmsscc2julian(raw_time);
    nRec = length(time);
    nc{'time'}(1:nRec) = floor(time);                                 
    nc{'time2'}(1:nRec) = (time-floor(time)).*(24*3600*1000); 
    clear raw_time time
else
    error('Times on raw data files should be the same')
end

% Define sample number
sample = [1:nSamples]';
nc{'sample'}(:) = sample;
clear sample

% Define burst number
burst = [1:nRec]';
nc{'burst'}(:) = burst;
clear burst

% Write lat/long to NetCDF file
nc{'lat'}(1)=nc.latitude(:);
nc{'lon'}(1)=nc.longitude(:);

% Write raw data to the raw data netCDF file
theVars = var(nc);
for v = 1:length(theVars),
    if strcmp(ncnames(theVars{v}),'press') || ...
       strcmp(ncnames(theVars{v}),'strk') || ...
       strcmp(ncnames(theVars{v}),'vel')
            theVar = char(ncnames(theVars{v}));
            eval(['data = burstData.',theVar,'.data;'])
            for rec = 1:nRec
                switch theVar
                    case 'press'
                        % reshape the press array so that
                        % rows is time and columns is sample
                        theVars{v}(rec, 1:nSamples) = data(:,rec)';
                        burstData.press.data = [];
                    
                    case 'strk'
                        % reshape the strk array so that
                        % rows is time, columns is sample, and
                        % the 3rd dimension is number of beams (4)
                        temp = reshape(data(:,rec), nBeams, length(data(:,rec))/nBeams)';
                        theVars{v}(rec, 1:nSamples, 1:nBeams) = temp;
                        clear temp
                        burstData.strk.data = [];
                        
                    case 'vel'
                        % reshape the vel array so that
                        % rows is time, columns is sample, and
                        % the 3rd dimension is number of beams (4)
                        % times number of bins used for waves (3),
                        % ergo, nBeamBins=12
                        temp = reshape(data(:,rec), nBeamBins, length(data)/nBeamBins)';
                        theVars{v}(rec, 1:nSamples, 1:nBeamBins) = temp;
                        clear temp burstData
                        burstData.vel.data = [];
                        
                    end
            end
            disp(['Finished writing ',theVar])
    end
end

% Calculate min/max values
calc_min_max_vals(nc);

% Update the history attribute
history = nc.history(:);
history_new = ['Pressure and velocity timeseries converted ',...
               'to NetCDF by adcpWvs2nc:ncwrite_adcpWvs_raw.m V ',...
               version,' on ',datestr(now,0),'; ',history];
nc.history = ncchar(history_new);

nc=close(nc);

disp(['Finished writing pressure and velocities. ',...
       num2str(toc/60), ' minutes elapsed'])
   
return


% ------------------------ Subfunctions --------------------------------- %
function [jd] = YYYYMMDDhhmmsscc2julian(YYYYMMDDhhmmsscc)
    
YYYY=floor(YYYYMMDDhhmmsscc/1000000000000);
MM=floor((YYYYMMDDhhmmsscc-YYYY.*1000000000000)./10000000000 );
DD=floor((YYYYMMDDhhmmsscc-(YYYY.*1000000000000)-(MM.*10000000000))./100000000);
hh=floor((YYYYMMDDhhmmsscc-(YYYY.*1000000000000)-(MM.*10000000000)-(DD.*100000000))./1000000);
mm=floor((YYYYMMDDhhmmsscc-(YYYY.*1000000000000)-(MM.*10000000000)-(DD.*100000000)-(hh.*1000000))./10000);
sec=floor((YYYYMMDDhhmmsscc-(YYYY.*1000000000000)-(MM.*10000000000)-(DD.*100000000)-(hh.*1000000)-(mm.*10000))./100);
hsec=floor((YYYYMMDDhhmmsscc-(YYYY.*1000000000000)-(MM.*10000000000)-(DD.*100000000)-(hh.*1000000)-(mm.*10000)-(sec.*100))); %1/100ths second

jd = julian([YYYY, MM, DD, hh, mm, sec+(hsec/100)]);

return

function calc_min_max_vals(nc)

theVars = var(nc);
for v = 1:length(theVars),
    if strcmp(ncnames(theVars{v}),'press') || ...
            strcmp(ncnames(theVars{v}),'strk') || ...
            strcmp(ncnames(theVars{v}),'vel')
        % calculate min/max values
        data = theVars{v}(:);
        theVars{v}.minimum = min(data);
        theVars{v}.maximum = max(data);
        clear data
    end
end

return
