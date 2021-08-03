function [specData]=read_adcpWvs_spec
% read_adcpWvs_spec.m  A function to load a series of ASCII files, which 
%                      are output from RD Instrument's WavesMon
%                      software, and output the directional and non-directional
%                      wave energy spectra.  
%
%   usage [specData]=read_adcpWvs_spec
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
%               
% Written by Charlene Sullivan
% USGS Woods Hole Science Center
% Woods Hole, MA 02543
% csullivan@usgs.gov

% C. Sullivan   10/26/05,   version 1.1
% Provide the user additional feedback regarding code execution. Get the
% direction at the start of the first direction slice in order to create a
% vector that is direction.  This is useful for plotting the directional
% spectra.
% C. Sullivan   06/09/05,   version 1.0
% This function assumes you are in a directory with WavesMon-generated wave
% data and that the files D-,P-,S-,VSpec*.txt exist in the directory. Bad
% data (0) is replaced with NaN. Don't do time conversion to julian
% days in here. Do that conversion when writing the data to NetCDF.


version = '1.1';

disp(' ')
disp(['Reading wave energy spectra'])

% Spectra types
specType = ['D', 'P', 'S', 'V'];


for s = 1:length(specType)
    %get list of files for the spectra type
    D = dir([specType(s),'Spec*.txt']);             
    nFiles = length(D);
       
    %loop through the files and load the data.  Also replace
    %all values of 0 (WavesMon bad data indicator for spectra
    %data) with NaN
    for n = 1:nFiles
        filename = D(n).name;
        
        if strcmp(filename(6),'0')
            filetime = str2num(filename(6:end-4)) + 200000000000;
        else
            filetime = str2num(filename(6:end-4)) + 190000000000;
        end
        
        switch specType(s)
            case 'D'
                if n == 1
                   %get the direction at which the first direction slice
                   %begins. this is determined by WavesMon and is the same
                   %throughout the deployment, but can vary between
                   %individual deployments
                   fid = fopen(filename,'r');
                   junk1 = fgetl(fid); junk2 = fgetl(fid); junk3 = fgetl(fid);
                   junk4 = fgetl(fid); junk5 = fgetl(fid); info = fgetl(fid);
                   clear junk*
                   fclose(fid);
                   specData.Dspec.firstDirSlice = ...
                       sscanf(info', '%*s %*s %*s %*s %*s %*s %*s %d %*s');                  
                end
                data = textread(filename,'%n','headerlines',6);
                data( data == 0 ) = nan;
                specData.Dspec.time(n) = filetime;
                specData.Dspec.data(:,n)=data;
            otherwise
                data = textread(filename,'%n','headerlines',3);
                data( data == 0 ) = nan;
                eval(['specData.',specType(s),'spec.time(n) = filetime;'])
                eval(['specData.',specType(s),'spec.data(:,n)=data;'])
        end
    end
    disp(['Finished reading wave energy spectra from ',specType(s),...
          'Spec*.txt. '])
end

return