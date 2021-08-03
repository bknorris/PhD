function nccreate_adcpWvs(metaFile,outFileRoot)
% nccreate_adcpWvs.m  A function to create empty netCDF files that will 
%                     store RD Instruments ADCP wave data.   
%
%    usage:  nccreate_adcpWvs(metaFile,outFileRoot);
%
%        where:  metaFile    - a string specifying the ascii file in which
%                              metadata is defined, in single quotes
%                              excluding the .txt file extension
%                outFileRoot - a string specifying the name given to the
%                              NetCDF output files, in single quotes
%                              excluding the NetCDF file extension .nc
%
% Written by Charlene Sullivan
% USGS Woods Hole Science Center
% Woods Hole, MA 02543
% csullivan@usgs.gov

% C. Sullivan   03/28/06,   version 1.2
% Add EPIC keys for the variables: maximum wave height (new EPIC key 4064),
% peak wave period (use existing EPIC key 4063), and peak wave direction
% (use existing EPIC key 4062). Don't use EPIC keys for spectra variables,
% because EPIC isn't suited for the spectral domain.  Perhaps CF conventions
% are better suited?
% C. Sullivan   10/25/05,   version 1.1
% Provide user additional feedback regarding code execution. File extension on
% metadata file no longer required in the input metaFile.  Changed DATA_TYPE
% attribute description to be consistant with the documentation.  Changed
% EPIC code and units on the variable lon to 502 and degree_east for
% consistancy w/ longitude as specified in the metafile where west is
% negative. Add ADCP bin size attributes.  Add selected bins for dspec and
% vspec to attributes. Define direction variable for directional spectra. 
% Re-define frequency variable.
% C. Sullivan   06/09/05,   version 1.0
% Now defining both a raw data and processed data NetCDF files.  Including
% all informaton from the waves configuration file as metadata.  Dimensions
% depth (=1), lat (=1), and lon (=1), beam (= 4), and beambin (=12) are
% hardwired.


version = '1.2';

ncr = netcdf([outFileRoot,'r-cal.nc'],'clobber');
ncp = netcdf([outFileRoot,'p-cal.nc'],'clobber');

% Gather and write NetCDF metadata
[userMeta, userMetaDefs] = textread([metaFile,'.txt'],'%s %63c','commentstyle','shell');
userMetaDefs = cellstr(userMetaDefs);

[wvmnMeta, wvmnMetaDefs] = read_WavesMon_config; 

write_adcpWvs_meta(ncr, userMeta, userMetaDefs, wvmnMeta, wvmnMetaDefs)
ncr.DATA_TYPE = ncchar('ADCP pressure and velocity timeseries');
ncr.VAR_DESC = ncchar('press:vel:strk');

write_adcpWvs_meta(ncp,userMeta, userMetaDefs, wvmnMeta, wvmnMetaDefs)
ncp.DATA_TYPE = ncchar('ADCP processed wave parameters and spectra');
ncp.VAR_DESC = ncchar('Hs:Tp:Dp:Hmax:Tm:dspec:pspec:sspec:vspec ');

% Define NetCDF dimensions
disp(['Defining NetCDF dimensions in ',outFileRoot,'r-cal.nc and ',...
       outFileRoot,'p-cal.nc'])
define_adcpWvs_dims(ncr);
define_adcpWvs_dims(ncp);

% Define NetCDF variables
disp(['Defining NetCDF variables in ',outFileRoot,'r-cal.nc and ',...
       outFileRoot,'p-cal.nc'])
disp(' ')
define_adcpWvs_vars(ncr);
define_adcpWvs_vars(ncp);

endef(ncr);
ncr=close(ncr);
endef(ncp);
ncp=close(ncp);

return

% --------- Subfunction: Gather WavesMon configuration data ------------- %
function [wvmnMeta,wvmnMetaDefs]=read_WavesMon_config;

    disp(' ')
    disp('Reading WavesMon configuration data ...')
    configFile = ls('*Wvs.cfg');
    ind = 1;
    fid = fopen(configFile,'r');
    while 1 
        tline = fgetl(fid);
        eqpos = strfind(tline,'=');
        if ~isempty(eqpos)
            wvmnMeta{ind,1} = tline(1:eqpos-1);
            wvmnMetaDefs{ind,1} = tline(eqpos+1:end);
            %convert whitespace in wvmnMeta to underscores
            ws = find(isspace(wvmnMeta{ind})==1); 
            wvmnMeta{ind}(ws) = '_';
            %convert parenthesis in parameter names to underscores
            op = strfind(wvmnMeta{ind},'(');
            wvmnMeta{ind}(op) = '_';
            cp = strfind(wvmnMeta{ind},')');
            wvmnMeta{ind}(cp) = '_';
            %convert commas in parameter names to underscores
            c = strfind(wvmnMeta{ind},',');
            wvmnMeta{ind}(c) = '_';
            %make sure no parameter name is longer than 63 characters
            %or Matlab complains
            if length(wvmnMeta{ind})>63
                wvmnMeta{ind} = wvmnMeta{ind}(1:63);
            end
            ind = ind + 1;
        end
        if ~ischar(tline), break, end
    end
    fclose(fid);
     
return

% --------- Subfunction: Write NetCDF metadata -------------------------- %
function write_adcpWvs_meta(nc, userMeta, userMetaDefs, wvmnMeta, wvmnMetaDefs)

    for a=1:length(userMeta)
        theField = userMeta{a};
        theFieldDef = userMetaDefs{a};
        if str2num(theFieldDef)
            nc.(theField) = str2num(theFieldDef);
        else
            nc.(theField) = theFieldDef;
        end
    end
    
    for a=1:length(wvmnMeta)
        theField = wvmnMeta{a};
        theFieldDef = wvmnMetaDefs{a};
        if str2num(theFieldDef) & ~isequal(theField,'VSpecBins')
            nc.WavesMonCfg.(theField) = str2num(theFieldDef);
        else
            nc.WavesMonCfg.(theField) = theFieldDef;
        end
    end
        
return
    
% --------- Subfunction: Define NetCDF dimensions ----------------------- %
function define_adcpWvs_dims(nc)

    %common dimensions for both raw and processed NetCDF files
    nc('time') = 0;
    nc('lat') = 1;
    nc('lon') = 1;
    
    if strcmp(nc.DATA_TYPE(:),'ADCP processed wave parameters and spectra')
        %dimensions specific to processed NetCDF file
        nc('frequency') = nc.WavesMonCfg.NFreqBins(:);
        nc('direction')  = nc.WavesMonCfg.NDir(:);
    elseif strcmp(nc.DATA_TYPE(:),'ADCP pressure and velocity timeseries')
        %dimensions specific to raw NetCDF file
        nc('sample') = nc.WavesMonCfg.FFTLen(:);
        nc('beam') = 4; %ADCP has 4 beams
        nc('beambin') = 4*3; %WavesMon uses 3 bins for each beam
    end

return

% --------- Subfunction: Define NetCDF variables ------------------------ %
function define_adcpWvs_vars(nc)
 
    %coordinate variables are the same for both the raw and
    %processed NetCDF files
    nc{'time'} = nclong('time'); 
    nc{'time'}.FORTRAN_format = ncchar('F10.2');
    nc{'time'}.units = ncchar('True Julian Day');
    nc{'time'}.type = ncchar('EVEN');
    nc{'time'}.epic_code = nclong(624);

    nc{'time2'} = nclong('time');
    nc{'time2'}.FORTRAN_format = ncchar('F10.2');
    nc{'time2'}.units = ncchar('msec since 0:00 GMT');
    nc{'time2'}.type = ncchar('EVEN');
    nc{'time2'}.epic_code = nclong(624);
    
    nc{'burst'} = nclong('time');
    nc{'burst'}.FORTRAN_format = ncchar('F10.2');
    nc{'burst'}.units = ncchar('counts');
    nc{'burst'}.type = ncchar('EVEN');

    nc{'lat'} = ncfloat('lat');
    nc{'lat'}.FORTRAN_format = ncchar('F10.4');
    nc{'lat'}.units = ncchar('degree_north');
    nc{'lat'}.type = ncchar('EVEN');
    nc{'lat'}.epic_code = nclong(500);

    nc{'lon'} = ncfloat('lon');
    nc{'lon'}.FORTRAN_format = ncchar('F10.4');
    nc{'lon'}.units = ncchar('degree_east');
    nc{'lon'}.type = ncchar('EVEN');
    nc{'lon'}.epic_code = nclong(502);

    
    if strcmp(nc.DATA_TYPE(:),'ADCP processed wave parameters and spectra')
        %Record variables for ONLY the processed data NetCDF file.
        %EPIC variables
        nc{'wh_4061'} = ncdouble('time','lat','lon');
        nc{'wh_4061'}.long_name=ncchar('Significant Wave Height (m)');
        nc{'wh_4061'}.generic_name=ncchar('wave_height');
        nc{'wh_4061'}.units=ncchar('m');
        nc{'wh_4061'}.epic_code=nclong(4061);
        nc{'wh_4061'}.FORTRAN_format=ncchar('F10.2');
        nc{'wh_4061'}.FillValue_ = ncfloat(1.0000000409184788E35);
        nc{'wh_4061'}.minimum = ncfloat(0);
        nc{'wh_4061'}.maximum = ncfloat(0);

        nc{'wp_4060'} = ncdouble('time','lat','lon');
        nc{'wp_4060'}.long_name=ncchar('Mean Wave Period (s)');
        nc{'wp_4060'}.generic_name=ncchar('wave_period');
        nc{'wp_4060'}.units=ncchar('s');
        nc{'wp_4060'}.epic_code=nclong(4060);
        nc{'wp_4060'}.FORTRAN_format=ncchar('F10.2');
        nc{'wp_4060'}.FillValue_ = ncfloat(1.0000000409184788E35);
        nc{'wp_4060'}.minimum = ncfloat(0);
        nc{'wp_4060'}.maximum = ncfloat(0);

        nc{'mwh_4064'} = ncdouble('time','lat','lon');
        nc{'mwh_4064'}.long_name=ncchar('Maximum Wave Height (m)');
        nc{'mwh_4064'}.generic_name=ncchar('wave_height');
        nc{'mwh_4064'}.units=ncchar('m');
        nc{'mwh_4064'}.epic_code=nclong(4064);
        nc{'mwh_4064'}.FORTRAN_format=ncchar('F10.2');
        nc{'mwh_4064'}.FillValue_ = ncfloat(1.0000000409184788E35);
        nc{'mwh_4064'}.minimum = ncfloat(0);
        nc{'mwh_4064'}.maximum = ncfloat(0);
        
        nc{'hght_18'} = ncdouble('time','lat','lon');
        nc{'hght_18'}.long_name=ncchar('Height of the Sea Surface (m)');
        nc{'hght_18'}.generic_name=ncchar('height');
        nc{'hght_18'}.units=ncchar('m');
        nc{'hght_18'}.epic_code=nclong(18);
        nc{'hght_18'}.FORTRAN_format=ncchar('F10.2');
        nc{'hght_18'}.NOTE=ncchar('height of sea surface relative to sensor');
        nc{'hght_18'}.FillValue_ = ncfloat(1.0000000409184788E35);
        nc{'hght_18'}.minimum = ncfloat(0);
        nc{'hght_18'}.maximum = ncfloat(0);

        nc{'wp_peak'} = ncdouble('time','lat','lon');
        nc{'wp_peak'}.long_name=ncchar('Peak Wave Period (s)');
        nc{'wp_peak'}.generic_name=ncchar('wave_period');
        nc{'wp_peak'}.units=ncchar('s');
        nc{'wp_peak'}.epic_code=nclong(4063);
        nc{'wp_peak'}.FORTRAN_format=ncchar('F10.2');
        nc{'wp_peak'}.FillValue_ = ncfloat(1.0000000409184788E35);
        nc{'wp_peak'}.minimum = ncfloat(0);
        nc{'wp_peak'}.maximum = ncfloat(0);

        nc{'wvdir'} = ncdouble('time','lat','lon');
        nc{'wvdir'}.long_name=ncchar('Peak Wave Direction (degrees true North)');
        nc{'wvdir'}.generic_name=ncchar('wave_dir');
        nc{'wvdir'}.units=ncchar('degrees T');
        nc{'wvdir'}.epic_code=nclong(4062);
        nc{'wvdir'}.FORTRAN_format=ncchar('F10.2');
        nc{'wvdir'}.NOTE=ncchar('Direction FROM which waves are propagating, measured clockwise from true north');
        nc{'wvdir'}.FillValue_ = ncfloat(1.0000000409184788E35);
        nc{'wvdir'}.minimum = ncfloat(0);
        nc{'wvdir'}.maximum = ncfloat(0);
        
        %non-EPIC variables
        nc{'frequency'}=ncdouble('frequency');
        nc{'frequency'}.name=ncchar('freq');
        nc{'frequency'}.long_name=ncchar('Frequency (Hz)');
        nc{'frequency'}.generic_name=ncchar('frequency');
        nc{'frequency'}.units=ncchar('Hz');
        nc{'frequency'}.type = ncchar('EVEN');
        nc{'frequency'}.FORTRAN_format=ncchar('F10.2');
        nc{'frequency'}.NOTE=ncchar('frequency at the center of each frequency band');
        nc{'frequency'}.minimum = ncfloat(0);
        nc{'frequency'}.maximum = ncfloat(0);
        
        nc{'direction'}=ncdouble('direction');
        nc{'direction'}.name=ncchar('dir');
        nc{'direction'}.long_name=ncchar('Direction (degrees T)');
        nc{'direction'}.generic_name=ncchar('direction');
        nc{'direction'}.units=ncchar('degrees T');
        nc{'direction'}.type = ncchar('EVEN');
        nc{'direction'}.FORTRAN_format=ncchar('F10.2');
        nc{'direction'}.NOTE=ncchar('direction at center of each direction slice');
        nc{'direction'}.minimum = ncfloat(0);
        nc{'direction'}.maximum = ncfloat(0);
        
        nc{'dspec'}=ncshort('time','frequency','direction','lat','lon');
        nc{'dspec'}.name=ncchar('dspec');
        nc{'dspec'}.long_name=ncchar('Directional Wave Energy Spectrum (mm^2/Hz/degree)');
        nc{'dspec'}.generic_name=ncchar('directional spectrum');
        nc{'dspec'}.units=ncchar('mm^2/Hz/degree');
        nc{'dspec'}.DspecBins=nc.WavesMonCfg.DirSpecBins(:);
        nc{'dspec'}.bin_size=nc.ADCPBinSize(:);
        nc{'dspec'}.FORTRAN_format=ncchar('F10.2');
        nc{'dspec'}.FillValue_ = ncshort(-32768);
        nc{'dspec'}.minimum = ncfloat(0);
        nc{'dspec'}.maximum = ncfloat(0);

        nc{'pspec'}=ncshort('time','frequency','lat','lon');
        nc{'pspec'}.name=ncchar('pspec');
        nc{'pspec'}.long_name=ncchar('Pressure-derived Non-directional Wave Height Spectrum (mm/sqrt(Hz))');
        nc{'pspec'}.generic_name=ncchar('pressure spectrum');
        nc{'pspec'}.units=ncchar('mm/sqrt(Hz)');
        nc{'pspec'}.FORTRAN_format=ncchar('F10.2');
        nc{'pspec'}.FillValue_ = ncshort(-32768);
        nc{'pspec'}.minimum = ncfloat(0);
        nc{'pspec'}.maximum = ncfloat(0);
        
        nc{'sspec'}=ncshort('time','frequency','lat','lon');
        nc{'sspec'}.name=ncchar('sspec');
        nc{'sspec'}.long_name=ncchar('Surface-derived Non-directional Wave Height Spectrum (mm/sqrt(Hz))');
        nc{'sspec'}.generic_name=ncchar('surface spectrum');
        nc{'sspec'}.units=ncchar('mm/sqrt(Hz)');
        nc{'sspec'}.FORTRAN_format=ncchar('F10.2');
        nc{'sspec'}.FillValue_ = ncshort(-32768);
        nc{'sspec'}.minimum = ncfloat(0);
        nc{'sspec'}.maximum = ncfloat(0);

        nc{'vspec'}=ncshort('time','frequency','lat','lon');
        nc{'vspec'}.name=ncchar('vspec');
        nc{'vspec'}.long_name=ncchar('Velocity-derived Non-directional Wave Height Spectrum (mm/sqrt(Hz))');
        nc{'vspec'}.generic_name=ncchar('velocity spectrum');
        nc{'vspec'}.units=ncchar('mm/sqrt(Hz)');
        nc{'vspec'}.VspecBins=nc.WavesMonCfg.VSpecBins(:);
        nc{'vspec'}.bin_size=nc.ADCPBinSize(:);
        nc{'vspec'}.FORTRAN_format=ncchar('F10.2');
        nc{'vspec'}.FillValue_ = ncshort(-32768);
        nc{'vspec'}.minimum = ncfloat(0);
        nc{'vspec'}.maximum = ncfloat(0);
        
    elseif strcmp(nc.DATA_TYPE(:),'ADCP pressure and velocity timeseries')
        %Record variables for ONLY the raw data NetCDF file.  There are no
        %EPIC-compatible variables for these quantities.
        nc{'sample'} = nclong('sample');
        nc{'sample'}.FORTRAN_format = ncchar('F10.2');
        nc{'sample'}.units = ncchar('counts');
        nc{'sample'}.type = ncchar('EVEN');
        
        nc{'press'}=ncshort('time','sample','lat','lon');
        nc{'press'}.name=ncchar('press');
        nc{'press'}.long_name=ncchar('Pressure Sensor Derived Depth (mm)');
        nc{'press'}.generic_name=ncchar('pressure time series');
        nc{'press'}.units=ncchar('mm');
        nc{'press'}.FORTRAN_format=ncchar('F10.2');
        nc{'press'}.FillValue_ = ncshort(-32768);
        nc{'press'}.minimum = ncfloat(0);
        nc{'press'}.maximum = ncfloat(0);
        
        nc{'strk'}=ncshort('time','sample','beam','lat','lon');
        nc{'strk'}.name=ncchar('strk');
        nc{'strk'}.long_name=ncchar('Along-Beam Surface Track (mm)');
        nc{'strk'}.generic_name=ncchar('surface track time series');
        nc{'strk'}.units=ncchar('mm');
        nc{'strk'}.FORTRAN_format=ncchar('F10.2');
        nc{'strk'}.FillValue_ = ncshort(-32768);
        nc{'strk'}.minimum = ncfloat(0);
        nc{'strk'}.maximum = ncfloat(0);
        
        nc{'vel'}=ncshort('time','sample','beambin','lat','lon');
        nc{'vel'}.name=ncchar('vspec');
        nc{'vel'}.long_name=ncchar('Along-Beam Velocity (mm/s)');
        nc{'vel'}.generic_name=ncchar('velocity time series');
        nc{'vel'}.units=ncchar('mm/s');
        nc{'vel'}.bin_size=nc.ADCPBinSize(:);
        nc{'vel'}.FORTRAN_format=ncchar('F10.2');
        nc{'vel'}.FillValue_ = ncshort(-32768);
        nc{'vel'}.minimum = ncfloat(0);
        nc{'vel'}.maximum = ncfloat(0);
        
    end
    
return
        

 
 