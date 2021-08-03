function Data = cropvp(Data,start,stop)
%function to crop a Vectrino Profiler file to the designated start and stop
%times provided by the user.
%
% Inputs: the Data structure of any VecPro file. Start and Stop must be in
% the matlab DATENUM format, i.e. datenum(yyyy,mm,dd,hh,mm,ss)
% Outputs: the Data structure with timeseries fields cropped to the start
% and stop times.
%
% This script was written by Benjamin K Norris, 2015
% University of Waikato, New Zealand
disp('Running Vectrino Quality Control')

if nargin < 3
	disp('Cannot crop timeseries without start and stop times')
    help(mfilename);
	return
end

disp(['Deployment Start time: ' datestr(start,'dd-mm-yyyy HH:MM:SS')])
disp(['Deployment Stop time: ' datestr(stop,'dd-mm-yyyy HH:MM:SS')])
disp('Cropping data to designated start/stop times')
gmt2ict = (datenum(0,0,0,1,0,0)*7); %ICT = GMT+7; VecPros record in GMT
if ~isfield(Data,'Time')
    fn = 1:22;
    Data.Profiles_HostTimeMatlab = Data.Profiles_HostTimeMatlab(1:end-1)+gmt2ict;
    ind = find(Data.Profiles_HostTimeMatlab >= start & Data.Profiles_HostTimeMatlab <= stop);
else
    fn = 1:23;
    Data.Time = Data.Time(1:end-1)+gmt2ict;
    ind = find(Data.Time >= start & Data.Time <= stop);
end
dfn = fieldnames(Data);
for i = fn
    disp(['Cropping field: ' dfn{i}])
    Data.(dfn{i}) = Data.(dfn{i})(ind,:);   
end
