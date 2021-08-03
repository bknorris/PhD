function aqdpWvs2nc(metaFile,outFileRoot);
% aqdpWvs2nc.m  A driver M-file for post-processing Nortek Aquadopp wave
%               data.
%
%    usage:  adcpWvs2nc(metaFile,outFileRoot);
%
%        where:  metaFile    - a string specifying the ascii file in which
%                              metadata is defined, in single quotes
%                              excluding the .txt file extension
%                outFileRoot - a string specifying the name given to the
%                              NetCDF output files, in single quotes
%                              excluding the NetCDF file extension .nc
%
% Written by Charlene Sullivan
% USGS Woods Hole Field Center
% Woods Hole, MA 02543
% csullivan@usgs.gov
%
% Toolbox Functions:
%   get_meta_nortek.m
%   wad2puv.m
%   nccreate_aqdpWvs.m
%   ncwrite_aqdpWvs.m
%   ncwrite_aqdpWvs_raw.m
%
% Add-on Functions:
%   beam2enu.m
%   wds.m
%   hs.m
%   julian.m
%   gmean.m
%   gmin.m
%   gmax.m

% C. Sullivan   10/27/05,   version 1.2
% Remove isunix loop b/c dir works on both unix and windows. Provide user
% additional feedback regarding code execution.  Metadata file must be a
% .txt file.  Don't require the .txt file extension in the input 'metaFile'.
% C. Sullivan   06/07/05,   version 1.1
% I also need the name of the .whd file output from AquaPro.  This file
% has heading, pitch, and roll information used to convert velocities
% from BEAM to ENU coordinates for PUV analysis.
% C. Sullivan   06/02/05,   version 1.0
% This Aquadopp component of the Wave Data Processing System toolbox writes
% both processed data and raw data NetCDF files. All user-defined metadata
% and instrument information is included in these files. The raw data
% NetCDF file contains the raw pressures and velocities used in the PUV
% analysis. User interaction is required to define which raw pressures and
% velocities are good (ie: in the water).  These good pressures and
% velocities are transformed from beam to geographic coordinates (if
% necessary) and then run through PUV analysis.  The results are written to
% the processed data NetCDF file.  Therefore, out-of-water data is excluded
% from PUV analysis and the processed data NetCDF file.


version = '1.2';

more off

tic

% Check inputs
if ~ischar(metaFile) || ~ischar(outFileRoot)
    error('File names should be surrounded in single quotes');
end

% Check existence of metadata file and Aquapro v. 1.25 output
l = dir([metaFile,'.txt']);
if isempty(l)
    error(['The metafile ',metaFile,'.txt does not exist in this directory']);
end

wad = dir('*.wad');
whd = dir('*.whd');
hdr = dir('*.hdr');
if isempty(wad) || isempty(whd) || isempty(hdr)
    error(['AquaPro output files do not exist in this directory'])
elseif length(wad) > 1 length(whd) > 1 || length(hdr) > 1
    error(['Too many .wad or .hdr files exist in this directory'])
else
    wadFile = wad.name;
    whdFile = whd.name;
    hdrFile = hdr.name;
end

% Gather user-defined and instrument metadata
[userMeta, aqdpMeta] = get_meta_nortek(metaFile, hdrFile);

% Perform PUV analysis
[burstData, puvData, aqdpMeta] = wad2puv(wadFile, whdFile, aqdpMeta, userMeta);

% Create and define your netcdf files
nccreate_aqdpWvs(userMeta, aqdpMeta, outFileRoot);

% Write data from puv analysis to processed data netcdf file
ncwrite_aqdpWvs(puvData, outFileRoot);

% Write raw pressures and velocities to raw data netcdf file
ncwrite_aqdpWvs_raw(burstData, outFileRoot);
