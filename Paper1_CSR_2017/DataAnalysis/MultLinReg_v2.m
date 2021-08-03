%Run the multiple linear regression function stepwise for the 1m^2 and
%0.2m^2 quadrat subsections. Look for significance in the model from
%various sources of turbulence: wave height, water depth, vegetation
%density, number of stems, average diameter, average velocity, z/h, and z/hc.

clear
%%%%Generate Table of hc values%%%%
path = 'd:\Projects\Mekong_W2015\DataAnalysis\Paper1\Environment\';
load([path 'VegetationData'])
hc = zeros(54,4);
hc(1,:) = vegdat.four.hc(1,1);hc(2,:) = vegdat.four.hc(1,2);hc(3,:) = vegdat.four.hc(1,3);
hc(4,:) = vegdat.four.hc(2,1);hc(5,:) = vegdat.four.hc(2,2);hc(6,:) = vegdat.four.hc(2,3);
hc(7,:) = vegdat.four.hc(3,1);hc(8,:) = vegdat.four.hc(3,2);hc(9,:) = vegdat.four.hc(3,3);
hc(10,:) = vegdat.four.hc(4,1);hc(11,:) = vegdat.four.hc(4,2);hc(12,:) = vegdat.four.hc(4,3);
hc(13,:) = vegdat.four.hc(5,1);hc(14,:) = vegdat.four.hc(5,2);hc(15,:) = vegdat.four.hc(5,3);
hc(16,:) = vegdat.four.hc(16,1);hc(17,:) = vegdat.four.hc(16,2);hc(18,:) = vegdat.four.hc(16,3);
hc(19,:) = vegdat.four.hc(17,1);hc(20,:) = vegdat.four.hc(17,2);hc(21,:) = vegdat.four.hc(17,3);
hc(22,:) = vegdat.four.hc(18,1);hc(23,:) = vegdat.four.hc(18,2);hc(24,:) = vegdat.four.hc(18,3);
hc(25,:) = vegdat.five.Hc(1,1);hc(26,:) = vegdat.five.Hc(1,2);hc(27,:) = vegdat.five.Hc(1,3);
hc(28,:) = vegdat.five.Hc(2,1);hc(29,:) = vegdat.five.Hc(2,2);hc(30,:) = vegdat.five.Hc(2,3);
hc(31,:) = vegdat.five.Hc(3,1);hc(32,:) = vegdat.five.Hc(3,2);hc(33,:) = vegdat.five.Hc(3,3);
hc(34,:) = vegdat.five.Hc(4,1);hc(35,:) = vegdat.five.Hc(4,2);hc(36,:) = vegdat.five.Hc(4,3);
hc(37,:) = vegdat.five.Hc(4,1);hc(38,:) = vegdat.five.Hc(4,2);hc(39,:) = vegdat.five.Hc(4,3);
hc(40,:) = vegdat.five.Hc(5,1);hc(41,:) = vegdat.five.Hc(5,2);hc(42,:) = vegdat.five.Hc(5,3);
hc(43,:) = vegdat.five.Hc(6,1);hc(44,:) = vegdat.five.Hc(6,2);hc(45,:) = vegdat.five.Hc(6,3);
hc(46,:) = vegdat.five.Hc(7,1);hc(47,:) = vegdat.five.Hc(7,2);hc(48,:) = vegdat.five.Hc(7,3);
hc(49,:) = vegdat.five.Hc(8,1);hc(50,:) = vegdat.five.Hc(8,2);hc(51,:) = vegdat.five.Hc(8,3);
hc(52,:) = vegdat.five.Hc(9,1);hc(53,:) = vegdat.five.Hc(9,2);hc(54,:) = vegdat.five.Hc(9,3);
%%%%
clear vegdat
load([path 'VelWvTKEdata'])
E = [veldat.four.E; veldat.five.E];
std = [veldat.four.Estd; veldat.five.Estd];
H = [veldat.four.depth; veldat.five.depth];
Hs = [veldat.four.wrms; veldat.five.wrms].*sqrt(2);
Avg = [veldat.four.avg; veldat.five.avg];

vpdir = 'd:\Projects\Mekong_W2015\DataAnalysis\DataReports\Paper1_QuadAnalysis\ThirdAttempt\';
files = dir([vpdir '*local20cm*.mat']);files = {files.name};
order = [2 4 3 1]; %load files in order
n = zeros(54,4);D = zeros(54,4);a = zeros(54,4);
for i = 1:length(files)
    load([vpdir files{order(i)}])
    stage = regexp(files{order(i)},'.+_(.*).mat','tokens');
    fname = ['twe' char(stage{:})];
    n(:,i) = vegdat.n;
    D(:,i) = vegdat.meanD;
    a(:,i) = vegdat.a;
end
z = vegdat.Zshore;

%%%Stepwise Functions%%%
%low tide
% zhc = z./hc(:,1);zhc(isinf(zhc)) = 0;
% data = [E(:,1) Hs(:,1) H(:,1) Avg(:,1) z./H(:,1) zhc a(:,1) n(:,1) D(:,1)];
% data(isnan(data)) = 0;
% X = [ones(length(E),1) data(:,2:end)];
% Y = data(:,1);
% stepwise(X,Y)

%mid-low tide
% Eps = E(:,2);Eps(Eps > 1E-3) = 0;
% zhc = z./hc(:,2);zhc(isinf(zhc)) = 0;
% data = [Eps Hs(:,2) H(:,2) Avg(:,2) z./H(:,2) zhc a(:,2) n(:,2) D(:,2)];
% data(isnan(data)) = 0;
% X = [ones(length(E),1) data(:,2:end)];
% Y = data(:,1);
% stepwise(X,Y)
% 
%mid-high tide
% zhc = z./hc(:,3);zhc(isinf(zhc)) = 0;
% data = [E(:,3) Hs(:,3) H(:,3) Avg(:,3) z./H(:,3) zhc a(:,3) n(:,3) D(:,3)];
% data(isnan(data)) = 0;
% X = [ones(length(E),1) data(:,2:end)];
% Y = data(:,1);
% stepwise(X,Y)
% 
%high tide
zhc = z./hc(:,4);zhc(isinf(zhc)) = 0;
data = [E(:,4) Hs(:,4) H(:,4) Avg(:,4) z./H(:,4) zhc a(:,4) n(:,4) D(:,4)];
data(isnan(data)) = 0;
X = [ones(length(E),1) data(:,2:end)];
Y = data(:,1);
stepwise(X,Y)

% %%%1m%%%
% load([vpdir 'Vegdat_1m.mat'])
% n = vegdat.n;
% D = vegdat.meanD;
% a = vegdat.a;
% avg = nanmean(Avg,2);
% h = nanmean(H,2);
% hs = nanmean(Hs,2);
% eps = nanmean(E,2);
% % 
% % %%%%
% zhc = z./hc(:,1);zhc(isinf(zhc)) = 0;
% data = [eps hs h avg z./h zhc a n D];
% X = [ones(length(E),1) data(:,2:end)];
% Y = data(:,1);
% stepwise(X,Y)
