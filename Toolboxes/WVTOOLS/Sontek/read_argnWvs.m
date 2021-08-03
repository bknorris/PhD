function [argnData]=read_argnWvs(datFile)
% read_aqdpWvs.m  A function to read wave data from a Sontek Argonaut.  
%
%   usage:  [argnData]=read_argnWvs;
%
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
%
% Written by Charlene Sullivan
% USGS Woods Hole Field Center
% Woods Hole, MA 02543
% csullivan@usgs.gov

% C. Sullivan   11/02/05,   version 1.1
% Provide user additional feedback regarding code execution.
% C. Sullivan   06/01/05,   version 1.0
% Output the data in a structure.  The data output in the structure will be
% in the same units as it is in the .dat file.  All unit conversions and
% calculation of the wave energy spectra will take place when the data is
% written to NetCDF.


version = '1.1';

argnData.datFile = datFile;

disp(['Reading statistical wave parameters from ',argnData.datFile])

% Load the data
data = load(argnData.datFile);
argnData.burst = [1:size(data,1)]';
argnData.YYYY = data(:,1);
argnData.MM = data(:,2);
argnData.DD = data(:,3);
argnData.hh = data(:,4);
argnData.mm = data(:,5);
argnData.ss = data(:,6);
argnData.Hs = data(:,47); 
argnData.Tp = data(:,48); 
argnData.ht = data(:,30);
argnData.ht_std = data(:,31);
argnData.amp = data(:,37:46);

return