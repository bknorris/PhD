clear
% vpdir = 'd:\Projects\Documents\Writing\DataReports\';
% fid = fopen([vpdir 'WaveStatsByDeployment.csv']);
% header = fgetl(fid);header = textscan(header,'%s','delimiter',',');
% pos = textscan(fid,'%s%s%s%n%n','delimiter',',');
% CMpos = cell2struct(pos,header{1},2);
% 
% save([vpdir 'CurrentMeterPos'],'CMpos','-v7.3')

vpdir = 'd:\Projects\Documents\Writing\DataReports\ThirdAttempt\';
files = dir([vpdir '*.csv']);files = {files.name};
for i = 1:length(files)
    fid = fopen([vpdir files{i}]);
    header = fgetl(fid);header = textscan(header,'%s','delimiter',',');
    pos = textscan(fid,'%s%s%n%n%n%n%n%n%n%n%n','delimiter',',');
    vegdat = cell2struct(pos,header{1},2);
    fname = regexprep(files{i},'.csv','');
    save([vpdir fname],'vegdat','-v7.3')
end

files = 'VPPositions.csv';
fid = fopen([vpdir files]);
header = fgetl(fid);header = textscan(header,'%s','delimiter',',');
pos = textscan(fid,'%s%s%n%n%n','delimiter',',');
VPpos = cell2struct(pos,header{1},2);
save([vpdir 'VP_Positions'],'VPpos','-v7.3')
