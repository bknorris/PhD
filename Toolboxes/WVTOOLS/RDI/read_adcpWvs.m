function [logData]=read_adcpWvs
% read_adcpWvs.m  A function to load *LogData.000, which is output from RD
%                 Instrument's WavesMon software, and output the
%                 timeseries of wave parameters. 
%
%   usage:  [logData]=read_adcpWvs;
%
%       where:  logData - a structure with the following fields
%                   file    - name of the *_LogData.* file
%                   burst   - burst number
%                   YY      - 2-digit year
%                   MM      - month
%                   DD      - day
%                   hh      - hours
%                   mm      - minutes
%                   ss      - seconds
%                   cc      - 1/100ths seconds
%                   Hs      - significant wave height, meters
%                   Hm      - maximum wave height, meters
%                   Tp      - peak wave period, seconds
%                   Tm      - mean wave period, seconds
%                   Dp      - peak wave direction, degrees true north
%                   ht      - water depth from pressure sensor, millimeters
%
% Written by Charlene Sullivan
% USGS Woods Hole Science Center
% Woods Hole, MA 02543
% csullivan@usgs.gov

% C. Sullivan   10/26/05,   version 1.1
% Provide user additional feedback regarding code execution.
% C. Sullivan   06/09/05,   version 1.0
% This function assumes you are in a directory with WavesMon-generated wave
% data and that the file *_LogData.000 exists in the directory.  Bad data
% indicator (-1) is replaced with NaN. Don't do time conversion to julian
% days in here. Do that conversion when writing the data to NetCDF.


version = '1.1';

% Load the *_LogData.* file
logFile = dir('*LogData*');
logData.file = logFile.name;
disp(['Reading statistical wave parameters from ',logData.file])
data = csvread(logData.file);

% Burst numbers
logData.burst = data(:,1);

% Extract time
logData.YY = data(:,2)+2000;
logData.MM = data(:,3);
logData.DD = data(:,4);
logData.hh = data(:,5);
logData.mm = data(:,6);
logData.ss = data(:,7);
logData.cc = data(:,8);

% Extract wave parameters
logData.Hs = data(:,9);
logData.Tp = data(:,10);
logData.Dp = data(:,11);
logData.ht = data(:,12);
logData.Hm = data(:,13);
logData.Tm = data(:,14);

% Replace all values of -1 (WavesMon bad data indicator
% for data is the *LogData.000 file) with NaN
theVars = {'Hs','Tp','Dp','ht','Hm','Tm'};
for v = 1:length(theVars)
    eval(['logData.',theVars{v},'(logData.',theVars{v},'== -1) = NaN;'])
end

return