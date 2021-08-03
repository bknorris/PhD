% Plot coastline and bathy for Molokai: Coastal Reefs Project

% Load in the location data
clear

fid = fopen('HW_outline.txt');
str = fgetl(fid);
HWhdr = textscan(str,'%s','Delimiter',',');
fields = HWhdr{1}';

FRMT = '%n%n%n%n%n';
C = textscan(fid,FRMT,'Delimiter',',');

HW = cell2struct(C,fields,2);

fid = fopen('SMKK_bathy.txt');
str = fgetl(fid);
MKKhdr = textscan(str,'%s','Delimiter',',');
fields = MKKhdr{1};

MKK = textscan(fid,FRMT,'Delimiter',',');

MKKbathy = cell2struct(MKK,fields,2);