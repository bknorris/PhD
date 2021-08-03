%Cat V5108 files together to make an environmental figure. Save U V P and
%time.
clear
datdir = 'c:\Users\bkn5\Projects\Mekong_W2015\Data\Vector\FSS\';
savedir = 'c:\Users\bkn5\Projects\Mekong_W2015\DataAnalysis\Environmental\';
A1 = load([datdir 'VC101_070315.mat']);
A2 = load([datdir 'VC101_100315.mat']);
% A3 = load([datdir 'V5108_120315.mat']);
e1 = A1.ADV.datetime(end);t2 = A2.ADV.datetime(1);
% e2 = A2.ADV.datetime(end);t3 = A3.ADV.datetime(1);
fs = 32;
int = datenum(0,0,0,0,0,1/fs);
tgap1 = (e1:int:t2)';
dgap1 = NaN(length(tgap1),1);
e1 = A1.ADV.Sensor.Datetime(end);t2 = A2.ADV.Sensor.Datetime(1);
% e2 = A2.ADV.datetime(end);t3 = A3.ADV.datetime(1);
fs = 1;
int = datenum(0,0,0,0,0,1/fs);
tgap2 = (e1:int:t2)';
dgap2 = NaN(length(tgap2),1);
% tgap2 = (e2:int:t3)';
% dgap2 = NaN(length(tgap2),1);
DATA = struct();
DATA.Metadata = A1.ADV.Metadata;
Time = [A1.ADV.datetime; tgap1; A2.ADV.datetime];% tgap2; A3.ADV.datetime];
Pres = [A1.ADV.Pres; dgap1; A2.ADV.Pres];% dgap2; A3.ADV.Pres];
Time2 = [A1.ADV.Sensor.Datetime; tgap2; A1.ADV.Sensor.Datetime];
Temp = [A1.ADV.Sensor.Temp; dgap2; A1.ADV.Sensor.Temp];
U = [A1.ADV.U; dgap1; A2.ADV.U];% dgap2; A3.ADV.U];
V = [A1.ADV.V; dgap1; A2.ADV.V];% dgap2; A3.ADV.V];
clear A1 A2 %A3

%use cmgbridge to fill gaps
nlin = 1E5;
maxg = 6E6;
DATA.datetime = Time;
DATA.Sensor.Datetime = Time2;
DATA.Sensor.Temp = cmgbridge(Temp,nlin,maxg,maxg);
DATA.Pres = cmgbridge(Pres,nlin,maxg,maxg);
% DATA.U = cmgbridge(U,nlin,maxg,maxg);
% DATA.V = cmgbridge(V,nlin,maxg,maxg);
DATA.U = U;
DATA.V = V;

ADV = DATA;
clear DATA
save([savedir 'VC101_HTA'],'-v7.3')
disp('File Saved')