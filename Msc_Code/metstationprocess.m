%Load and process MET station weather data
clear
dirname = 'c:\Users\bkn5\Projects\Mekong_F2014\Data\Weather\';
fname = 'Cu Lao Dung-ST0-ST1-Weather-2014.csv';
fid = fopen([dirname fname]);
str = fgetl(fid);
header = textscan(str,'%s','Delimiter',',');
fields = header{1};
FRMT = '%s%n%n%n%n%n%n%s';
thedata = textscan(fid,FRMT,'delimiter',',');
Met = cell2struct(thedata,fields,2);

datetime = NaN(length(Met.DateTime),1);
%Convert time into serial date no
for i = 1:length(Met.DateTime)
    if isempty(Met.DateTime{i})
        continue
    else
    datetime(i) = datenum(Met.DateTime{i},'dd/mm/yyyy HH:MM');
    end
end
Met = rmfield(Met,'DateTime');
Met.Datetime = datetime;


