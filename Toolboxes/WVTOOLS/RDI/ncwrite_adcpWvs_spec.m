function ncwrite_adcpWvs_spec(specData, outFileRoot)
% ncwrite_adcpWvs_spec.m  A function to write the timeseries of directional
%                         and non-directional spectra, output from RD
%                         Instruments WavesMon, to a netCDF file.  
%
%   usage:  ncwrite_adcpWvs_spec(specData, outFileRoot);
%
%       where:  specData - a structure with the following fields
%                   dspec.data  - directional wave energy spectra, mm^2/(Hz)/deg
%                   dspec.time  - dspec time, YYYYMMDDhhmm   
%                   pspec.data  - pressure-derived non-directional wave energy
%                                 spectra, mm/sqrt(Hz) 
%                   pspec.time  - pspec time, YYYYMMDDhhmm 
%                   sspec.data  - surface-derived non-directional wave energy
%                                 spectra, mm/sqrt(Hz)
%                   sspec.time  - sspec time, YYYYMMDDhhmm 
%                   vspec.data  - velocity-derived non-directional wave energy
%                                 spectra, mm/sqrt(Hz)
%                   vspec.time  - vspec time
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

% C. Sullivan   03/29/06,   version 1.3
% Reverse the chronology on the history attribute so the most recent
% processing step is listed first. In calculating max and min attributes
% for the spectra, calculate the max and min over time for each frequency
% bin.
% C. Sullivan   03/08/06,   version 1.2
% The function no longer assumes the times of the spectra data include
% seconds. When aligning spectra data with wave parameters just use year,
% month, day, hour, and minute.
% C. Sullivan   10/26/05,   version 1.1
% Provide the user more feedback regarding code execution.  Define the
% variable direction for the directional spectra.  This direction is the
% direction at the center of each directional slice for the directional
% spectra.
% C. Sullivan   06/09/05,   version 1.0
% This function assumes times on the spectra data include the 30 seconds
% that the times in the *_LogData* files have. 


version = '1.3';

nc = netcdf([outFileRoot,'p-cal.nc'],'write');

% Get dimensions from NetCDF file
theRecDim = nc('time');
nRec = size(theRecDim);
nRec = nRec(1);
theFreqDim = nc('frequency');
nFreq = size(theFreqDim);
nFreq = nFreq(1);
dFreq = 1/nFreq; % Frequency bands are 1/nFreq Hz wide
theDirDim = nc('direction');
nDir = size(theDirDim);
nDir = nDir(1);
dDir = 360/nDir; % Direction slices are nDir degrees wide

% Define a vector that is frequency.  This vector consists of the frequency
% at the CENTER of each frequency band.
frequency = [0:dFreq:dFreq*(nFreq-1)]';
frequency = (frequency + (frequency+dFreq))/2;

% Define a vector that is direction for the directional spectra.  This
% vector consists of the direction in the CENTER of each direction slice.
firstDir = specData.Dspec.firstDirSlice;
direction = [firstDir : dDir : dDir*(nDir-1)+firstDir]';
direction = (direction + (direction+dDir-1))/2;
f = find(direction > 360);
direction(f) = direction(f) - 360;

% Get time from the NetCDF file
nc_time = gregorian(nc{'time'}(:) + nc{'time2'}(:)/3600/1000/24);
nc_time = julian(nc_time);

% Write spectra data to the NetCDF file
disp(' ')
disp(['Writing wave energy spectra to ',outFileRoot,'p-cal.nc'])
theVars = var(nc);
for v = 1:length(theVars),
    if strcmp(ncnames(theVars{v}),'frequency')  || ...
       strcmp(ncnames(theVars{v}),'direction')  || ...
       strcmp(ncnames(theVars{v}),'dspec') || ...
       strcmp(ncnames(theVars{v}),'pspec') || ...
       strcmp(ncnames(theVars{v}),'sspec') || ...
       strcmp(ncnames(theVars{v}),'vspec')
            if strcmp(ncnames(theVars{v}),'frequency')
                % write frequency to NetCDF file
                theVars{v}(1:nFreq) = frequency;
                theVars{v}.minimum = gmin(frequency);
                theVars{v}.maximum = gmax(frequency);
            elseif strcmp(ncnames(theVars{v}),'direction')
                % write direction to NetCDF file
                theVars{v}(1:nDir) = direction;
                theVars{v}.minimum = gmin(direction);
                theVars{v}.maximum = gmax(direction); 
            else
                if strcmp(ncnames(theVars{v}),'dspec')
                    % reshape the dspec data array so that
                    % rows is time, columns is frequency, and the
                    % 3rd dimension is direction. this time dimension
                    % differs from that in the NetCDF file, so I'm 
                    % initializing another array to which I will later align
                    % (in time) the dspec data.
                    data_new = nan(nRec, nFreq, nDir);
                    sz = size(specData.Dspec.data);
                    for i = 1:sz(2)
                        temp = reshape(specData.Dspec.data(:,i), nDir, nFreq)';
                        data(i,:,:) = temp;
                    end
                    spec_time = specData.Dspec.time';
                else
                    % reshape the p(s)(v)spec data arrays so that
                    % rows is time, and columns is frequency. this time dimension
                    % differs from that in the NetCDF file, so I'm 
                    % initializing another array to which I will later align
                    % (in time) the p(s)(v)spec data.
                    data_new = nan(nRec, nFreq);
                    spec = ncnames(theVars{v});
                    spec = upper(spec{1}(1));
                    eval(['data = specData.',spec,'spec.data'';'])
                    eval(['spec_time = specData.',spec,'spec.time'';']);
                end
                % convert times on the spectra from YYYYMMDDhhmm to
                % julian. then align the spectra data (in time) with the
                % time records in the NetCDF file. the time records in the
                % NetCDF file are based on times in the *LogData.000 file 
                % (the wave parameters file).
                spec_time = YYYYMMDDhhmm2julian(spec_time);
                for rec=1:nRec
                    srec1 = find(spec_time>=nc_time(rec),1,'first');
                    dist1 = spec_time(srec1) - nc_time(rec);
                    srec2 = find(spec_time<=nc_time(rec),1,'last');
                    dist2 = nc_time(rec) - spec_time(srec2);
                    if ~isempty(srec1) && ~isempty(srec2)
                        if srec1 == srec2
                            srec = srec1;
                        elseif dist1 > dist2
                            srec = srec2;
                        elseif dist2 > dist1
                            srec = srec1;
                        elseif dist1 == dist2
                            srec = srec1;
                        end
                    else
                        srec = [];
                    end
                    if ~isempty(srec)
                        if strcmp(ncnames(theVars{v}),'dspec')
                            data_new(rec, :, :) = data(srec, :, :);
                        else
                            data_new(rec, :) = data(srec, :);
                        end
                    end
                end
                
                clear data spec_time
                
                % write the spectra data to NetCDF file and replace any
                % nans in the data with the NetCDF fill value
                if strcmp(ncnames(theVars{v}),'dspec')
                    theFillVal = theVars{v}.FillValue_(:);
                    bads = find(isnan(data_new));
                    data_new(bads) = nan;
                    theVars{v}.minimum = gmin(data_new(:));
                    theVars{v}.maximum = gmax(data_new(:));
                    data_new(bads) = theFillVal;
                    theVars{v}(1:nRec, 1:nFreq, 1:nDir) = data_new;
                else
                    theFillVal = theVars{v}.FillValue_(:);
                    bads = find(isnan(data_new));
                    data_new(bads) = nan;
                    theVars{v}.minimum = gmin(data_new);
                    theVars{v}.maximum = gmax(data_new);
                    data_new(bads) = theFillVal;
                    theVars{v}(1:nRec, 1:nFreq) = data_new;
                end
                theVar = char(ncnames(theVars{v}));
                disp(['Finished writing ',theVar])
                clear data_new
            end
    end                
end                     

% Update the history attribute
history = nc.history(:);
history_new = ['Wave energy spectra converted to NetCDF by ',...
               'adcpWvs2nc:ncwrite_adcpWvs_spec V ',version,...
               ' on ',datestr(now,0),'; ',history];
nc.history = ncchar(history_new);

nc=close(nc);

disp(['Finished writing wave energy spectra. ',...
      num2str(toc/60),' minutes elapsed'])

return

% ------------------------ Subfunctions --------------------------------- %
function [jd] = YYYYMMDDhhmm2julian(YYYYMMDDhhmm);

YYYY=floor(YYYYMMDDhhmm/100000000);
MM=floor((YYYYMMDDhhmm-YYYY.*100000000)./1000000 );
DD=floor((YYYYMMDDhhmm-(YYYY.*100000000)-(MM.*1000000))./10000);
hh=floor((YYYYMMDDhhmm-(YYYY.*100000000)-(MM.*1000000)-(DD.*10000))./100);
mm=floor((YYYYMMDDhhmm-(YYYY.*100000000)-(MM.*1000000)-(DD.*10000)-(hh.*100)));

jd = julian([YYYY, MM, DD, hh, mm, zeros(size(YYYY))]);

return
