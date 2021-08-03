%Big loadnplot for all the CLAY runs

clear
filelist = dir('*PAR.mat');
files = {filelist.name};
dirname = 'D:\Projects\PARSensor\Clay';
%Do the load thing
names = {'blk1000';'blk100';'blk10';'red1000';'red100';'red10';'tan1000';'tan100';'tan10';...
    'wht1000';'wht100';'wht10'}';
data = {};
for i = 1:length(files)
    fname = fullfile(dirname,files{i});
    data{i} = load(fname);
end
par = cell2struct(data,names,2);
depth = 0:5:160; %step interval (in cm)

filelist = dir('*SLOBS.mat');
files = {filelist.name};
%Do the load thing
names = {'blk1000';'blk100';'blk10';'red1000';'red100';'red10';'tan1000';'tan100';'tan10';...
    'wht1000';'wht100';'wht10'}';
data = {};
for i = 1:length(files)
    fname = fullfile(dirname,files{i});
    data{i} = load(fname);
end
slobs = cell2struct(data,names,2);


