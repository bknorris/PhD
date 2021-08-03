function nccreate_aqdpWvs(userMeta, aqdpMeta, outFileRoot)
% nccreate_aqdpWvs.m  A function to create empty netCDF files that will
%                     store Nortek Aquadopp wave data.   
%
%    usage:  nccreate_aqdpWvs(userMeta, aqdpMeta, outFileRoot);
%
%        where: userMeta - a structure with user-defined metadata
%               aqdpMeta - a structure with Aquadopp setup parameters
%               outFileRoot - a string specifying the name given to the
%                             NetCDF output files, in single quotes
%                             excluding the NetCDF file extension .nc
%
% Written by Charlene Sullivan
% USGS Woods Hole Field Center
% Woods Hole, MA 02543
% csullivan@usgs.gov

% C. Sullivan   03/29/06,   version 1.3
% Add EPIC codes to the variables wp_peak and wvdir.
% C. Sullivan   10/28/05,   version 1.2
% Provide the user additional feedback regarding code execution. Changed
% DATA_TYPE attribute description to be consistent w/ the documentation.
% Changed EPIC code and units on the variable lon to 502 and degree_east
% for consistency w/ longitude as specified in the metadata file where west
% is negative.
% C. Sullivan   06/09/05,   version 1.1
% The variables 'pspec' and 'vspec' will now have units of mm/sqrt(Hz) for
% consistency with RDI and SonTek spectra. (the old units were m^2/Hz)
% C. Sullivan   06/02/05,   version 1.0
% PUV analysis now run as part of the toolbox. Instrument setup information
% is now included in the metadata.  Two NetCDF output files are created;
% one file has the raw pressure and velocity data that was used in PUV
% analysis; the other has the processed data from the PUV analysis.


version = '1.2';

ncr = netcdf([outFileRoot,'r-cal.nc'],'clobber');
ncp = netcdf([outFileRoot,'p-cal.nc'],'clobber');

write_aqdpWvs_meta(ncr, userMeta, aqdpMeta)
ncr.DATA_TYPE = ncchar('Aquadopp pressure and velocity timeseries');
ncr.VAR_DESC = ncchar('u:v:w:p');

write_aqdpWvs_meta(ncp,userMeta, aqdpMeta)
ncp.DATA_TYPE = ncchar('Aquadopp processed wave parameters and spectra');
ncp.VAR_DESC = ncchar('Hs:Tp:Dp:Pspec:Vspec');

% Define NetCDF dimensions
disp(' ')
disp(['Defining NetCDF dimensions in ',outFileRoot,'r-cal.nc and ',...
       outFileRoot,'p-cal.nc'])
define_aqdpWvs_dims(ncr);
define_aqdpWvs_dims(ncp);

% Define NetCDF variables
disp(['Defining NetCDF variables in ',outFileRoot,'r-cal.nc and ',...
       outFileRoot,'p-cal.nc'])
disp(' ')
define_aqdpWvs_vars(ncr);
define_aqdpWvs_vars(ncp);

endef(ncr);
ncr=close(ncr);
endef(ncp);
ncp=close(ncp);

return

% --------- Subfunction: Write NetCDF metadata -------------------------- %
function write_aqdpWvs_meta(nc, userMeta, aqdpMeta)

    theAtts = fieldnames(userMeta);
    
    for a=1:length(theAtts)
        eval(['theDef = userMeta.',theAtts{a},';'])
        if ischar(theDef)
            eval(['nc.',theAtts{a},' = ncchar(''',theDef,''');']);
        else
            eval(['nc.',theAtts{a},' = ncfloat([theDef]);']);
        end
    end
    
    theAtts = fieldnames(aqdpMeta);
    
    for a=1:length(theAtts)
        eval(['theDef = aqdpMeta.',theAtts{a},';'])
        if ischar(theDef)
            eval(['nc.',theAtts{a},' = ncchar(''',theDef,''');']);
        else
            eval(['nc.',theAtts{a},' = ncfloat([theDef]);']);
        end
    end
        
return
    
% --------- Subfunction: Define NetCDF dimensions ----------------------- %
function define_aqdpWvs_dims(nc)

    %common dimensions for both raw and processed NetCDF files
    nc('time') = 0;
    nc('lat') = 1;
    nc('lon') = 1;
    
    if strcmp(nc.DATA_TYPE(:),'Aquadopp processed wave parameters and spectra')
        %dimensions specific to processed NetCDF file
        nc('frequency') = nc.nf_out(:); %this is a result from the PUV analysis
    elseif strcmp(nc.DATA_TYPE(:),'Aquadopp pressure and velocity timeseries')
        %dimensions specific to raw NetCDF file
        nc('sample') = nc.Wave__Number_of_samples(:);
    end

return

% --------- Subfunction: Define NetCDF variables ------------------------ %
function define_aqdpWvs_vars(nc)
 
    %coordinate variables are the same for both the raw and
    %processed NetCDF files
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

    
    if strcmp(nc.DATA_TYPE(:),'Aquadopp processed wave parameters and spectra')
        %Record variables for ONLY the processed data NetCDF file.
        %EPIC variables
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

        nc{'wh_4061'} = ncdouble('time','lat','lon');
        nc{'wh_4061'}.long_name=ncchar('Significant Wave Height (m)');
        nc{'wh_4061'}.generic_name=ncchar('Hs');
        nc{'wh_4061'}.units=ncchar('m');
        nc{'wh_4061'}.epic_code=nclong(4061);
        nc{'wh_4061'}.FORTRAN_format=ncchar('F10.2');
        nc{'wh_4061'}.FillValue_ = 1.0000000409184788E35;
        nc{'wh_4061'}.minimum = ncfloat(0);
        nc{'wh_4061'}.maximum = ncfloat(0);

        nc{'hght_18'} = ncdouble('time','lat','lon');
        nc{'hght_18'}.long_name=ncchar('Height of the Sea Surface (m)');
        nc{'hght_18'}.generic_name=ncchar('height');
        nc{'hght_18'}.units=ncchar('m');
        nc{'hght_18'}.epic_code=nclong(18);
        nc{'hght_18'}.FORTRAN_format=ncchar('F10.2');
        nc{'hght_18'}.NOTE=ncchar('height of sea surface relative to sensor');
        nc{'hght_18'}.FillValue_ = 1.0000000409184788E35;
        nc{'hght_18'}.minimum = ncfloat(0);
        nc{'hght_18'}.maximum = ncfloat(0);

        %non-EPIC variables
        nc{'frequency'}=ncdouble('frequency');
        nc{'frequency'}.name=ncchar('freq');
        nc{'frequency'}.long_name=ncchar('Frequency (Hz)');
        nc{'frequency'}.generic_name=ncchar('frequency');
        nc{'frequency'}.units=ncchar('Hz');
        nc{'frequency'}.type=ncchar('EVEN');
        nc{'frequency'}.NOTE=ncchar('frequency at the center of each frequency band');
        nc{'frequency'}.FORTRAN_format=ncchar('F10.2');
        nc{'frequency'}.minimum = ncfloat(0);
        nc{'frequency'}.maximum = ncfloat(0);
        
        nc{'dfreq'}=ncdouble('frequency');
        nc{'dfreq'}.name=ncchar('freq');
        nc{'dfreq'}.long_name=ncchar('Frequency Band Width (Hz)');
        nc{'dfreq'}.generic_name=ncchar('freq band width');
        nc{'dfreq'}.units=ncchar('Hz');
        nc{'dfreq'}.FORTRAN_format=ncchar('F10.2');
        nc{'dfreq'}.minimum = ncfloat(0);
        nc{'dfreq'}.maximum = ncfloat(0);
        nc{'dfreq'}.NOTE = ncchar('Bandwidth of each frequency band');
        
        nc{'wp_peak'} = ncdouble('time','lat','lon');
        nc{'wp_peak'}.long_name=ncchar('Peak Wave Period (s)');
        nc{'wp_peak'}.generic_name=ncchar('Tp');
        nc{'wp_peak'}.units=ncchar('s');
        nc{'wp_peak'}.epic_code=nclong(4063);
        nc{'wp_peak'}.FORTRAN_format=ncchar('F10.2');
        nc{'wp_peak'}.FillValue_ = 1.0000000409184788E35;
        nc{'wp_peak'}.minimum = ncfloat(0);
        nc{'wp_peak'}.maximum = ncfloat(0);

        nc{'wvdir'} = ncdouble('time','lat','lon');
        nc{'wvdir'}.long_name=ncchar('Peak Wave Direction (degrees true North)');
        nc{'wvdir'}.generic_name=ncchar('Dp');
        nc{'wvdir'}.units=ncchar('degrees');
        nc{'wvdir'}.epic_code=nclong(4062);
        nc{'wvdir'}.FORTRAN_format=ncchar('F10.2');
        nc{'wvdir'}.NOTE=ncchar('Direction FROM which waves are propagating');
        nc{'wvdir'}.FillValue_ = 1.0000000409184788E35;
        nc{'wvdir'}.minimum = ncfloat(0);
        nc{'wvdir'}.maximum = ncfloat(0);

        nc{'spread'} = ncdouble('time','lat','lon');
        nc{'spread'}.long_name=ncchar('Peak Spreading');
        nc{'spread'}.generic_name=ncchar('spreading');
        nc{'spread'}.units=ncchar('degrees');
        nc{'spread'}.FORTRAN_format=ncchar('F10.2');
        nc{'spread'}.FillValue_ = 1.0000000409184788E35;
        nc{'spread'}.minimum = ncfloat(0);
        nc{'spread'}.maximum = ncfloat(0);
        
        nc{'pspec'}=ncfloat('time','frequency','lat','lon');
        nc{'pspec'}.name=ncchar('pspec');
        nc{'pspec'}.long_name=ncchar('Pressure-derived Non-directional Wave Height Spectrum (mm/sqrt(Hz))');
        nc{'pspec'}.generic_name=ncchar('pressure spectrum');
        nc{'pspec'}.units=ncchar('mm/sqrt(Hz)');
        nc{'pspec'}.FORTRAN_format=ncchar('F10.2');
        nc{'pspec'}.FillValue_ = 1.0000000409184788E35;
        nc{'pspec'}.minimum = ncfloat(0);
        nc{'pspec'}.maximum = ncfloat(0);
        
        nc{'vspec'}=ncfloat('time','frequency','lat','lon');
        nc{'vspec'}.name=ncchar('vspec');
        nc{'vspec'}.long_name=ncchar('Velocity-derived Non-directional Wave Height Spectrum (mm/sqrt(Hz))');
        nc{'vspec'}.generic_name=ncchar('velocity spectrum');
        nc{'vspec'}.units=ncchar('mm/sqrt(Hz)');
        nc{'vspec'}.FORTRAN_format=ncchar('F10.2');
        nc{'vspec'}.FillValue_ = 1.0000000409184788E35;
        nc{'vspec'}.minimum = ncfloat(0);
        nc{'vspec'}.maximum = ncfloat(0);

    elseif strcmp(nc.DATA_TYPE(:),'Aquadopp pressure and velocity timeseries')
        %Record variables for ONLY the raw data NetCDF file.  There are no
        %EPIC-compatible variables for these quantities.
        nc{'time'} = nclong('time','sample'); 
        nc{'time'}.FORTRAN_format = ncchar('F10.2');
        nc{'time'}.units = ncchar('True Julian Day');
        nc{'time'}.type = ncchar('EVEN');
        nc{'time'}.epic_code = nclong(624);

        nc{'time2'} = nclong('time','sample');
        nc{'time2'}.FORTRAN_format = ncchar('F10.2');
        nc{'time2'}.units = ncchar('msec since 0:00 GMT');
        nc{'time2'}.type = ncchar('EVEN');
        nc{'time2'}.epic_code = nclong(624);
        
        nc{'sample'} = nclong('sample');
        nc{'sample'}.FORTRAN_format = ncchar('F10.2');
        nc{'sample'}.units = ncchar('counts');
        nc{'sample'}.type = ncchar('EVEN');
        
        nc{'hght_18'} = ncdouble('time','sample','lat','lon');
        nc{'hght_18'}.long_name=ncchar('Height of the Sea Surface (m)');
        nc{'hght_18'}.generic_name=ncchar('height');
        nc{'hght_18'}.units=ncchar('m');
        nc{'hght_18'}.epic_code=nclong(18);
        nc{'hght_18'}.FORTRAN_format=ncchar('F10.2');
        nc{'hght_18'}.NOTE=ncchar('height of sea surface relative to sensor');
        nc{'hght_18'}.FillValue_ = 1.0000000409184788E35;
        nc{'hght_18'}.minimum = ncfloat(0);
        nc{'hght_18'}.maximum = ncfloat(0);
        
        nc{'u_1205'}=ncfloat('time','sample','lat','lon');
        nc{'u_1205'}.name=ncchar('u');
        nc{'u_1205'}.long_name=ncchar('Eastward Velocity (cm/s)');
        nc{'u_1205'}.generic_name=ncchar('u');
        nc{'u_1205'}.units=ncchar('cm/s');
        nc{'u_1205'}.epic_code=nclong(1205);
        nc{'u_1205'}.FORTRAN_format=ncchar('');
        nc{'u_1205'}.FillValue_ = 1.0000000409184788E35;
        nc{'u_1205'}.minimum = ncfloat(0);
        nc{'u_1205'}.maximum = ncfloat(0);
        
        nc{'v_1206'}=ncfloat('time','sample','lat','lon');
        nc{'v_1206'}.name=ncchar('v');
        nc{'v_1206'}.long_name=ncchar('Northward Velocity (cm/s)');
        nc{'v_1206'}.generic_name=ncchar('v');
        nc{'v_1206'}.units=ncchar('cm/s');
        nc{'v_1206'}.epic_code=nclong(1206);
        nc{'v_1206'}.FORTRAN_format=ncchar('');
        nc{'v_1206'}.FillValue_ = 1.0000000409184788E35;
        nc{'v_1206'}.minimum = ncfloat(0);
        nc{'v_1206'}.maximum = ncfloat(0);
        
        nc{'w_1204'}=ncfloat('time','sample','lat','lon');
        nc{'w_1204'}.name=ncchar('w');
        nc{'w_1204'}.long_name=ncchar('Vertical Velocity (cm/s)');
        nc{'w_1204'}.generic_name=ncchar('w');
        nc{'w_1204'}.units=ncchar('cm/s');
        nc{'w_1204'}.epic_code=nclong(1204);
        nc{'w_1204'}.FORTRAN_format=ncchar('');
        nc{'w_1204'}.FillValue_ = 1.0000000409184788E35;
        nc{'w_1204'}.minimum = ncfloat(0);
        nc{'w_1204'}.maximum = ncfloat(0);
        
        nc{'amp'}=ncfloat('time','sample','lat','lon');
        nc{'amp'}.units=ncchar('counts');
        nc{'amp'}.long_name=ncchar('Beam amplitude');
        nc{'amp'}.FillValue_ = 1.0000000409184788E35;
    end
    
return