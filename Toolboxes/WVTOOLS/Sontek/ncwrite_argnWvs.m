function ncwrite_argnWvs(argnData, outFileRoot)
% ncwrite_adcpWvs.m  A function to write wave data from a Sontek Argonaut.  
%
%   usage:  ncwrite_aqdpWvs(argnData, outFileRoot);
%       where:  argnData - a structure with the following fields
%                   datFile - name of the .dat file to which ViewArgonaut
%                             output data
%                   burst   - burst number 
%                   YYYY    - 4-digit year
%                   MM      - month
%                   DD      - day
%                   hh      - hours
%                   mm      - minutes
%                   ss      - seconds
%                   Hs      - significant wave height, cm
%                   Tp      - peak wave period, s
%                   ht      - water depth from pressure sensor, dBar
%                   ht_std  - standard deviation water depth, dBar
%                   amp     - wave amplitude, mm
%              outFileRoot - a string specifying the name given to the
%                            NetCDF output files, in single quotes
%                            excluding the NetCDF file extension .nc
%
% Written by Charlene Sullivan
% USGS Woods Hole Field Center
% Woods Hole, MA 02543
% csullivan@usgs.gov
%
% Dependencies:
%   julian.m
%   gregorian.m
%   gmin.m
%   gmax.m

% C.Sullivan    03/30/06,   version 1.3
% Reverse the chronology on the history attribute so the most recent
% processing step is listed first.  For max and min attributes on the
% variable pspec, calculate the maxs and mins over time for each frequency
% band.
% C.Sullivan    11/02/05,   version 1.2
% Provide the user additional feedback regarding code execution.
% C. Sullivan   06/09/05,   version 1.1
% Users is interactively asked to view a plot to determine which bursts are
% good.  Only good bursts are written to NetCDF.
% C. Sullivan   06/01/05,   version 1.0
% Unit conversions and the calculation of wave spectra take place in this
% m-file.  Data is written to the processed data NetCDF file.  No raw data
% NetCDF file is produced as the 1 Hz pressure timeseries from which wave
% parameters are calculated is not output by ViewArgonaut.
 

version = '1.2';

% ---- Define good data to write to NetCDF ------------------------------ %
%display a figure with two subplots.  the first subplot has mean pressure
%and the standard deviation of pressure.  the second subplot contains the
%significant wave height, and the peak wave period.  the user is asked to
%view these plots (zooming in and out) and determine the first and last
%good bursts of data.  all data between the first and last
%good data points will be written to NetCDF.
F=figure;
set(F,'position',[200 100 800 800])
subplot(411)
plot(argnData.burst, argnData.ht)
ylabel({'Pressure','(dBar)'})
subplot(412)
plot(argnData.burst, argnData.ht_std)
ylabel({'Pressure Standard','Deviation (dBar)'})
subplot(413)
plot(argnData.burst, argnData.Hs)
ylabel({'Significant','Wave Height (cm)'})
subplot(414)
plot(argnData.burst, argnData.Tp)
ylabel({'Peak Wave','Period (s)'})
xlabel('Index')

disp(' ')
disp('Please view the figure (feel free to zoom!) and look for indices which might be bad')
disp('Bad indices might include out-of-water indices, for example')
answer = input(['Would you like to write all indices (indices 1:',...
              num2str(length(argnData.burst)),') to NetCDF? y/n :  '],'s');
if ~strcmp(answer,'Y') & ~strcmp(answer,'y') & ...
   ~strcmp(answer,'N') & ~strcmp(answer,'n')
        disp(['Valid answers are either yes ("y") or no ("n")'])
        answer = input(['Would you like to write all indices (indices 1:',...
                       num2str(length(argnData.burst)),') to NetCDF? y/n :  '],'s');
end
if strcmp(answer,'Y') || strcmp(answer,'y')
    first_good = 1;
    last_good = length(argnData.burst);
    disp(['Including indices 1 through ',num2str(length(argnData.burst)),' in PUV analysis'])
elseif strcmp(answer,'N') || strcmp(answer,'n')
    first_good = input('Please view the figure and enter the first good indice:  ');
    if isempty(first_good)
        disp('Valid answers are numeric only')
        first_good = input('Please view the figure and enter the first good indice:  ');
    end
    last_good = input('Please view the figure and enter the last good indice:  ');
    if isempty(last_good)
        disp('Valid answers are numeric only')
        last_good = input('Please view the figure and enter the last good indice:  ');
    end
    disp(['Writing indices ',num2str(first_good),' to ',num2str(last_good), ...
          ' to NetCDF'])
end
goodBursts = [first_good:last_good];
nGoodBursts = length(goodBursts);

theVars = fieldnames(argnData);
for v = 1:length(theVars)
    if ~strcmp(theVars{v},'datFile')
        if strcmp(theVars{v},'amp')
            argnData.amp = argnData.amp(goodBursts,:);
        else
            eval(['argnData.',theVars{v},' = argnData.',theVars{v},'(goodBursts);'])
        end
    end
end

% ---- Calculate wave energy density ------------------------------------ %
%convert 'amplitude at 10 period bands' into energy density in mm/sqrt(Hz)
%according to P. Work's equations.
%wave period bands are 2 seconds wide.  the first band begins at 2 seconds.
%the last bin is for periods > 20 seconds.  
T_low = [2:2:20]'; 
T_high = [4:2:22]'; 
T_high(end) = 22; %  <--- ??? larger period defining the last bin
dF = 1./T_low - 1./T_high;

totEnergy = 0.5 .* sum((argnData.amp.^2)); %total energy in each band
Hmo = 4 .* sqrt(totEnergy); % zero-moment wave height in each band

for i=1:length(dF)
    E(:,i) = (0.5*(argnData.amp(:,i).^2)) ./ dF(i); %energy density
                                                    %(cm^2/Hz) in each band
end

argnData.pspec = sqrt(E).*10; %we call pressure-derived energy density
                              %'pspec' and the units we like are mm/sqrt(Hz)
                     
argnData.frequency = 1./T_low + dF/2; %we also like to have the frequency
                                      %at the center of each frequency band

% ---- Write data to to NetCDF ------------------------------------- %

disp(' ')
disp(['Writing statistical wave parameters and spectra to ',outFileRoot,'p-cal.nc'])

nc=netcdf([outFileRoot,'p-cal.nc'],'write');

% Convert time from gregorian to julian and get
% start_time and stop_time
time = julian([argnData.YYYY, argnData.MM, argnData.DD,...
               argnData.hh, argnData.mm, argnData.ss]);
start_time=datenum([gregorian(time(1))]);
stop_time=datenum([gregorian(time(end))]);

% Convert Hs from cm to m
argnData.Hs = argnData.Hs ./ 100;

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
nc{'burst'}(1:nRec) = argnData.burst;
nc{'lat'}(1)=nc.latitude(:);
nc{'lon'}(1)=nc.longitude(:);

% write record variables
ncvarnames = {'wh_4061','wp_peak','hght_18','hght_std','frequency','pspec'};
names = {'Hs','Tp','ht','ht_std','frequency','pspec'};
nNames = length(names);
for i = 1:nNames
    eval(['data = argnData.',names{i},';'])
    eval(['nco = nc{''',ncvarnames{i},'''};'])
    nco.maximum = gmax(data);
    nco.minimum = gmin(data);
    theFillVal = nco.FillValue_(:);
    bads = find(isnan(data));
    data(bads) = theFillVal;
    if strcmp(ncvarnames(i),'pspec')
        nco(1:nRec, 1:nFreq) = data; 
    elseif strcmp(ncvarnames(i),'frequency')
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
history_new = ['Wave energy spectra calculated and data converted ',...
               'to NetCDF by argnWvs2nc:ncwrite_argnWvs.m V ',version,...
               ' on ' datestr(now,0),'; ',history];
nc.history = history_new;

% Close netCDF file
nc=close(nc);

elapsed_time = toc;
disp(['Finished writing statistical wave parameters and spectra. ',...
      num2str(toc/60),' minutes elapsed.'])
  
return