function nccreate_argnWvs(userMeta, argnMeta, outFileRoot)
% nccreate_argnWvs.m  A function to create empty netCDF files that will
%                     store wave data from a Sontek Argonaut.   
%
%    usage:  nccreate_argnWvs(userMeta, argnMeta, outFileRoot);
%
%        where: userMeta - a structure with user-defined metadata
%               argnMeta - a structure with Argonaut setup parameters
%               outFileRoot -  a string specifying the name given to the
%                             NetCDF output files, in single quotes
%                             excluding the NetCDF file extension .nc
%
% Written by Charlene Sullivan
% USGS Woods Hole Field Center
% Woods Hole, MA 02543
% csullivan@usgs.gov

% C. Sullivan   03/30/06,   version 1.2
% Use EPIC keys 1239 and 4063 for the variables hght_std and wp_peak,
% respectively.
% C. Sullivan   11/02/05,   version 1.1
% Provide the user additional feedback regarding code execution. Changed
% DATA_TYPE attribute description to be consistent w/ the documentation.
% Changed EPIC code and units on the variable lon to 502 and degree_east
% for consistency w/ longitude as specified in the metadata file where west
% is negative.
% C. Sullivan   06/01/05,   version 1.0
% Instrument setup information is now included in the metadata.  Only the
% processed data NetCDF file is created; the file includes wave parameters 
% which were calculated from a 1 Hz pressure timeseries.  No raw data
% NetCDF file is created as the 1 Hz pressure timeseries s not output by
% ViewArgonaut.  It is assumed there are 10 frequency bands in the wave
% energy spectra.


version = '1.2';

ncp = netcdf([outFileRoot,'p-cal.nc'],'clobber');

write_argnWvs_meta(ncp, userMeta, argnMeta)
ncp.DATA_TYPE = ncchar('Argonaut processed wave parameters and spectra');
ncp.VAR_DESC = ncchar('Hs:Tp:Pspec');

% Define NetCDF dimensions
disp(' ')
disp(['Defining NetCDF dimensions in ',outFileRoot,'p-cal.nc.'])
define_argnWvs_dims(ncp);

% Define NetCDF variables
disp(['Defining NetCDF variables in ',outFileRoot,'p-cal.nc.'])
disp(' ')
define_argnWvs_vars(ncp);

endef(ncp);
ncp=close(ncp);

return

% --------- Subfunction: Write NetCDF metadata -------------------------- %
function write_argnWvs_meta(nc, userMeta, argnMeta)

    theAtts = fieldnames(userMeta);
    
    for a=1:length(theAtts)
        eval(['theDef = userMeta.',theAtts{a},';'])
        if ischar(theDef)
            eval(['nc.',theAtts{a},' = ncchar(''',theDef,''');']);
        else
            eval(['nc.',theAtts{a},' = ncfloat([theDef]);']);
        end
    end
    
    theAtts = fieldnames(argnMeta);
    
    for a=1:length(theAtts)
        eval(['theDef = argnMeta.',theAtts{a},';'])
        if ischar(theDef)
            eval(['nc.',theAtts{a},' = ncchar(''',theDef,''');']);
        else
            eval(['nc.',theAtts{a},' = ncfloat([theDef]);']);
        end
    end
        
return
    
% --------- Subfunction: Define NetCDF dimensions ----------------------- %
function define_argnWvs_dims(nc)

    nc('time') = 0;
    nc('lat') = 1;
    nc('lon') = 1;
    nc('frequency') = 10; % it is assumed there are 10 frequency bands

return

% --------- Subfunction: Define NetCDF variables ------------------------ %
function define_argnWvs_vars(nc)
 
    %coordinate variables
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
    
    %record variables -- EPIC variables
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

    %record variables -- non-EPIC variables
    nc{'hght_std'} = ncdouble('time','lat','lon');
    nc{'hght_std'}.long_name=ncchar('Standard Deviation of Height of the Sea Surface (m)');
    nc{'hght_std'}.generic_name=ncchar('height std');
    nc{'hght_std'}.units=ncchar('m');
    nc{'hght_std'}.epic_code=nclong(1239);
    nc{'hght_std'}.FORTRAN_format=ncchar('F10.2');
    nc{'hght_std'}.FillValue_ = 1.0000000409184788E35;
    nc{'hght_std'}.minimum = ncfloat(0);
    nc{'hght_std'}.maximum = ncfloat(0);
    
    nc{'frequency'}=ncdouble('frequency');
    nc{'frequency'}.name=ncchar('freq');
    nc{'frequency'}.long_name=ncchar('Frequency (Hz)');
    nc{'frequency'}.generic_name=ncchar('frequency');
    nc{'frequency'}.units=ncchar('Hz');
    nc{'frequency'}.NOTE=ncchar('frequency at the center of each frequency band');
    nc{'frequency'}.FORTRAN_format=ncchar('F10.2');
    nc{'frequency'}.minimum = ncfloat(0);
    nc{'frequency'}.maximum = ncfloat(0);

    nc{'wp_peak'} = ncdouble('time','lat','lon');
    nc{'wp_peak'}.long_name=ncchar('Peak Wave Period (s)');
    nc{'wp_peak'}.generic_name=ncchar('Tp');
    nc{'wp_peak'}.units=ncchar('s');
    nc{'wp_peak'}.epic_code=nclong(4063);
    nc{'wp_peak'}.FORTRAN_format=ncchar('F10.2');
    nc{'wp_peak'}.FillValue_ = 1.0000000409184788E35;
    nc{'wp_peak'}.minimum = ncfloat(0);
    nc{'wp_peak'}.maximum = ncfloat(0);

    nc{'pspec'}=ncfloat('time','frequency','lat','lon');
    nc{'pspec'}.name=ncchar('pspec');
    nc{'pspec'}.long_name=ncchar('Pressure-derived Non-directional Wave Height Spectrum (mm/sqrt(Hz))');
    nc{'pspec'}.generic_name=ncchar('velocity spectrum');
    nc{'pspec'}.units=ncchar('mm/sqrt(Hz)');
    nc{'pspec'}.FORTRAN_format=ncchar('F10.2');
    nc{'pspec'}.FillValue_ = 1.0000000409184788E35;
    nc{'pspec'}.minimum = ncfloat(0);
    nc{'pspec'}.maximum = ncfloat(0);

return