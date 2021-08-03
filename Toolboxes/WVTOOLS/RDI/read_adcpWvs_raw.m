function [burstData]=read_adcpWvs_raw
% read_adcpWvs_raw.m  A function to load Press*.txt, STrk*.txt, and
%                     Vel*.txt which is output from RD Instrument's
%                     WavesMon v. 2.01 software, and output the pressure,
%                     surface track, and velocity timeseries. 
%
%   usage:  [burstData]=read_adcpWvs_raw;
%
%       where:  burstData - an structure w/ the following fields:
%                   burstData.press.data - the pressure timeseries,
%                   burstData.press.time - press time, YYYYMMDDhhmmsscc
%                   burstData.strk.data - the range-to-surface track timeseries,
%                   burstData.strk.time - strk time, YYYYMMDDhhmmsscc
%                   burstData.vel.data - the velocity timeseries 
%                   burstData.vel.time - vel time, YYYYMMDDhhmmsscc
%
% Written by Charlene Sullivan
% USGS Woods Hole Science Center
% Woods Hole, MA 02543
% csullivan@usgs.gov

% C. Sullivan   10/26/06,   version 1.1
% Provide the user more feedback regarding code execution.
% C. Sullivan   06/09/05,   version 1.0
% This function assumes you are in a directory with WavesMon-generated
% wave data and that the files Press*.txt, STrk*.txt, and Vel*.txt exist in
% the directory.


version = '1.1';

disp(' ')
disp('Reading pressure and velocity timeseries')

% Raw data types
rawType = {'Press', 'Strk', 'Vel'};

for s = 1:length(rawType)
    
    %get list of files for the raw data type  
    D = dir([rawType{s},'*.txt']);             
    nFiles = length(D);
       
    %loop through the files and load the data.
    for n = 1:nFiles
        filename = D(n).name;
        
        switch rawType{s}
            case 'Press'
                %load pressure timeseries
                filetime = str2num(filename(6:end-4));
                fid = fopen(filename,'r');
                fgetl(fid); fgetl(fid); fgetl(fid);
                data = fscanf(fid,'%f');
                burstData.press.time(n) = filetime;
                burstData.press.data(:,n) = int16(data);
                fclose(fid);
                
            case 'Strk'
                %load range-to-surface track timeseries
                filetime = str2num(filename(5:end-4));
                fid=fopen(filename,'r');
                fgetl(fid); fgetl(fid); fgetl(fid); fgetl(fid);
                data = fscanf(fid,'%f %f %f %f');
                burstData.strk.time(n) = filetime;
                burstData.strk.data(:,n) = int16(data);
                fclose(fid);
            
            case 'Vel'
                %load velocity timeseries
                filetime = str2num(filename(4:end-4));
                fid=fopen(filename,'r');
                fgetl(fid); fgetl(fid); fgetl(fid);
                fgetl(fid); fgetl(fid); fgetl(fid);
                data = fscanf(fid,'%f %f %f %f %f %f %f %f %f %f %f %f');
                burstData.vel.time(n) = filetime;
                burstData.vel.data(:,n) = int16(data);
                fclose(fid);
        end
    end
    disp(['Finished reading ',rawType{s}])
   
end    

disp(['Finished reading pressures and velocities'])

return