function argnWvs2nc(metaFile,outFileRoot);
% argnWvs2nc.m  A driver m-file for post-processing wave data from a
%               Sontek Argonaut.  
%
%    usage:  argnWvs2nc(metaFile,outFileRoot);
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
%   get_meta_sontek.m
%   nccreate_argnWvs.m
%   read_argnWvs.m
%   ncwrite_argnWvs.m
%
% Add-on Functions:
%   julian.m
%   gregorian.m
%   gmin.m
%   gmax.m

% C.Sullivan    11/02/05,   version 1.1
% Provide user additional feedback regarding code execution. File extension
% on metadata file no longer required for the input 'metaFile', but the
% metadata file must be a text file w/ the .txt file extension. Remove
% isunix loop b/c dir works on both Unix and Windows.
% C. Sullivan   06/01/05,   version 1.0
% This Argonaut component of the Wave Data Processing System toolbox only
% writes a processed data NetCDF file, as the 1 hz pressure timeseries
% from which wave parameters are calculated, is not output by ViewArgonaut.


more off

version = '1.1';

tic

% Check that input and output file names are characters
if ~ischar(metaFile) || ~ischar(outFileRoot)
    error('Input and output file names should be surrounded in single quotes');
end

% Check existence of metadata file and ViewArgonaut output.
l = dir([metaFile,'.txt']);
if isempty(l)
    error(['The metafile ',metaFile,'.txt does not exist in this directory']);
end

ctl = dir('*.ctl');
if isempty(ctl)
    error(['ViewArgonaut output files do not exist in this directory'])
else
    ctlFile = ctl.name;
    datFile = [ctlFile(1:end-4),'.dat'];
end

% Gather user-defined and instrument metadata
[userMeta, argnMeta] = get_meta_sontek(metaFile, ctlFile);

% Create and define your netcdf file
nccreate_argnWvs(userMeta, argnMeta, outFileRoot);

% Load data
[argnData]=read_argnWvs(datFile);

% Write data to netcdf file
ncwrite_argnWvs(argnData, outFileRoot);
