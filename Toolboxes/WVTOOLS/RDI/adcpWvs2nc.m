function adcpWvs2nc(metaFile,outFileRoot)
% adcpWvs2nc.m  A driver M-file for post-processing RD Instruments
%               WavesMon wave data in Matlab.  
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
% USGS Woods Hole Science Center
% Woods Hole, MA 02543
% csullivan@usgs.gov
%
% Toolbox Functions:
%   nccreate_adcpWvs.m
%   read_adcpWvs.m
%   ncwrite_adcpWvs.m
%   read_adcWvs_spec.m
%   ncwrite_adcpWvs_spec.m
%   read_adcpWvs_raw.m
%   ncwrite_adcpWvs_raw.m
%
% Add-on Functions:
%   julian.m
%   gregorian.m
%   gmin.m
%   gmax.m

% C. Sullivan   10/25/05,   version 1.1
% Provide user additional feedback regarding code execution. File extension
% on metadata file no longer required for the input 'metaFile'.
% C. Sullivan   06/09/05,   version 1.0
% This function assumes the user ran RD Instruments WavesMon 
% software, output the WavesMon raw and processed data to text files, and
% created the file 'metaFile.dat' with important metadata. The function
% assumes you are in a directory with this data and the netcdf toolbox and
% the Wave Data Processing System toolbox exist on your Matlab path. Raw
% and processed data is written to individual netCDF files.  An integer
% fill value (-32768) is used for the spectra data and the raw data.
% Clear statements interspersed throughout mfiles prevent out of memory
% errors.


more off

version = '1.1';

tic

% Check inputs
if ~ischar(metaFile) || ~ischar(outFileRoot)
    error('File names should be surrounded in single quotes');
end

% Check existence of metadata file
l=ls([metaFile,'.txt']);
if isempty(l)
    error(['The metafile ',metaFile,'.txt does not exist in this directory']);
end

% Check WavesMon output
if isempty(dir('*_LogData.*'))
    error(['WavesMon output does not exist in the directory ',pwd])
elseif length(dir('*_LogData.*')) > 1
    error('Only one *_LogData.* file is permitted by this toolbox')
elseif isempty(dir('*Spec*.txt'))
    error('WavesMon output must include spectra data as text files')
elseif isempty(dir('Strk*.txt'))
    error('WavesMon output must include raw Range to surface data')
elseif isempty(dir('Press*.txt'))
    error('WavesMon output must include raw pressure data')
elseif isempty(dir('Vel*.txt'))
    error('WavesMon output must include raw orbital velocity data') 
end
    
% Create and define your netcdf files
nccreate_adcpWvs(metaFile, outFileRoot);

% Load timeseries data
[logData] = read_adcpWvs;

% Write timeseries data to NetCDF
ncwrite_adcpWvs(logData, outFileRoot);

clear logData

% Load spectra data
[specData] = read_adcpWvs_spec;

% Write spectra data to NetCDF
ncwrite_adcpWvs_spec(specData, outFileRoot);

clear specData 

% Load raw data
[rawData] = read_adcpWvs_raw;

% Write raw data to NetCDF
ncwrite_adcpWvs_raw(rawData, outFileRoot);

clear rawData