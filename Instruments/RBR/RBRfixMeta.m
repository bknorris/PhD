%do a little organizing with the RBR logger data

dirc = 'C:\Users\bkn5\Projects\Mekong_F2014\Data\RBR\SoloP\Raw\';
% comment = 'Mekong 2014 Fine Scale Study 1 variant 2';
fileList = dir([dirc '*.mat']);
filelist = {fileList.name}';

%Solo P Gauges
for i = 1:length(filelist)
    disp(['Loading file ' filelist{i}])
    load([dirc filelist{i}])
    
    %assign missing fields in metadata
    RBR.Metadata.units = {'Pres: dBar';'SeaPres: dBar';'Depth: m'};
    RBR.Metadata.hab = [num2str(RBR.Metadata.hab) ' mm'];
    RBR.Metadata.Deployment_date = datestr(RBR.Metadata.deployment_date,'dd-mmm-yyyy HH:MM:SS');
    RBR.Metadata.Recovery_date = datestr(RBR.Metadata.recovery_date,'dd-mmm-yyyy HH:MM:SS');
    RBR.Metadata = rmfield(RBR.Metadata,'starttime');
    RBR.Metadata = rmfield(RBR.Metadata,'endtime');
    RBR.Metadata = rmfield(RBR.Metadata,'deployment_date');
    RBR.Metadata = rmfield(RBR.Metadata,'recovery_date');
    RBR.Metadata.nsamp = length(RBR.Datetime);

    disp('Saving File...')
    save([dirc filelist{i}],'RBR')
    
    clear RBR
end

%Duet dual wave & temperature gauge
% for i = 1:length(filelist)
%     disp(['Loading file ' filelist{i}])
%     load([dirc filelist{i}])
%     
%     %assign missing fields in metadata
%     RBR.Metadata.units = {'Temp: C';'Pres: dBar';'SeaPres: dBar';'Depth: m'};
%     RBR.Metadata.hab = [num2str(RBR.Metadata.hab) ' mm'];
%     RBR.Metadata.Deployment_date = datestr(RBR.Metadata.deployment_date,'dd-mmm-yyyy HH:MM:SS');
%     RBR.Metadata.Recovery_date = datestr(RBR.Metadata.recovery_date,'dd-mmm-yyyy HH:MM:SS');
%     RBR.Metadata = rmfield(RBR.Metadata,'starttime');
%     RBR.Metadata = rmfield(RBR.Metadata,'endtime');
%     RBR.Metadata = rmfield(RBR.Metadata,'deployment_date');
%     RBR.Metadata = rmfield(RBR.Metadata,'recovery_date');
%     RBR.Metadata.nsamp = length(RBR.Datetime);
%     
%     disp('Saving File...')
%     save([dirc filelist{i}],'RBR')
%     
%     clear RBR
% end

% %Solo T Loggers
% for i = 1:length(filelist)
%     disp(['Loading file ' filelist{i}])
%     load([dirc filelist{i}])
%     assign missing fields in metadata
%     RBR.Metadata.units = 'Degrees C';
%     RBR.Metadata.Deployment_date = datestr(RBR.Metadata.deployment_date,'dd-mmm-yyyy HH:MM:SS');
%     RBR.Metadata.Recovery_date = datestr(RBR.Metadata.recovery_date,'dd-mmm-yyyy HH:MM:SS');
%     RBR.Metadata = rmfield(RBR.Metadata,'deployment_date');
%     RBR.Metadata = rmfield(RBR.Metadata,'recovery_date');
%     RBR.Metadata = rmfield(RBR.Metadata,'chanunits');
%     RBR.Metadata.nsamp = length(RBR.Datetime);
%     
%     disp('Saving File...')
%     save([dirc filelist{i}],'RBR')
%     
%     clear RBR
% end

%fix the 2014 CTDs
% RBR.Datetime = RBR.datetime;RBR = rmfield(RBR,'datetime');
% RBR.Yearday = RBR.yearday;RBR = rmfield(RBR,'yearday');
% RBR.Metadata.instname = RBR.name;
% RBR.Metadata.Serial = 18816;
% RBR.Metadata.data_cmt = 'Mekong 2014 Dense Pneumatophore Study NE Cu Lao Dung';
% RBR.Metadata.samprate = '6 Hz';
% RBR.Metadata.channelnames = {'Conductivity';'Temperature';'Pressure';'Depth';'Salinity';'Specific Conductivity';'Density Anomaly';'Speed of sound'};
% RBR.Metadata.Deployment_date = '02-10-2014 11:00:00';
% RBR.Metadata.Recovery_date = '03-10-2014 12:00:00';
% RBR.Metadata.units = RBR.chanunits;
% RBR.Metadata.nsamp = length(RBR.Datetime);
% RBR = rmfield(RBR,'chanunits');
% RBR = rmfield(RBR,'name');
% RBR = rmfield(RBR,'samprate');
% RBR = rmfield(RBR,'param');

%fix the 2014 Solo T loggers
% for i = 1:length(filelist)
%     disp(['Loading file ' filelist{i}])
%     load([dirc filelist{i}])
%     
%     RBR.Datetime = RBR.datetime';RBR = rmfield(RBR,'datetime');
%     RBR.Yearday = RBR.yearday';RBR = rmfield(RBR,'yearday');
%     %assign missing fields in metadata
%     RBR.Metadata.inst_type = 'RBR Solo T Logger';
%     RBR.Metadata.data_cmt = comment;
%     RBR.Metadata.instname = RBR.name;
%     RBR.Metadata.serial = str2double(RBR.name(end-5:end));
%     RBR.Metadata.samprate = '4 Hz';
%     RBR.Metadata.channelnames = {'Temperature'};
%     RBR.Metadata.coefficients = [0.0035;-0.0003;0.0000;-0.0000];
%     RBR.Metadata.parameters = RBR.param;
%     RBR.Metadata.units = 'Degrees C';
%     RBR.Metadata.Deployment_date = datestr(RBR.Datetime(1),'dd-mmm-yyyy HH:MM:SS');
%     RBR.Metadata.Recovery_date = datestr(RBR.Datetime(end),'dd-mmm-yyyy HH:MM:SS');
%     RBR = rmfield(RBR,'name');
%     RBR = rmfield(RBR,'samprate');
%     RBR = rmfield(RBR,'param');
%     RBR = rmfield(RBR,'chanunits');
%     RBR.Metadata.nsamp = length(RBR.Datetime);
%     
%     disp('Saving File...')
%     save([dirc filelist{i}],'RBR')
%     
%     clear RBR
% end



disp('Finished!')