clear

%First, load the data from CatTurbData.m and CatVegGeometry.m
p1 = 'D:\Projects\Code\ImageProcessing\';
run([p1 'CatVegGeometry.m'])
p2 = 'D:\Projects\Code\DataAnalysis\Turbulence\';
run([p2 'CatTurbData_v2.m'])

%%%Update 09/06/2016: To plot phi over the entire area sampled (1m^2),
%%%uncomment the following section and comment the section at the bottom
%%%where phi appears.
% phi = zeros(24,3);
% phi(1:3,:) = repmat(vegdat.four.Phi(1,:),3,1);
% phi(4:6,:) = repmat(vegdat.four.Phi(2,:),3,1);
% phi(7:9,:) = repmat(vegdat.four.Phi(3,:),3,1);
% phi(10:12,:) = repmat(vegdat.four.Phi(4,:),3,1);
% phi(13:15,:) = repmat(vegdat.four.Phi(5,:),3,1);
% phi(16:18,:) = repmat(vegdat.four.Phi(16,:),3,1);
% phi(19:21,:) = repmat(vegdat.four.Phi(17,:),3,1);
% phi(22:24,:) = repmat(vegdat.four.Phi(18,:),3,1);
% phif = phi; clear phi
% phi = zeros(27,3);
% phi(1:3,:) = repmat(vegdat.five.Phi(1,:),3,1);
% phi(4:6,:) = repmat(vegdat.five.Phi(2,:),3,1);
% phi(7:9,:) = repmat(vegdat.five.Phi(3,:),3,1);
% phi(10:12,:) = repmat(vegdat.five.Phi(4,:),3,1);
% phi(13:15,:) = repmat(vegdat.five.Phi(5,:),3,1);
% phi(16:18,:) = repmat(vegdat.five.Phi(6,:),3,1);
% phi(19:21,:) = repmat(vegdat.five.Phi(7,:),3,1);
% phi(22:24,:) = repmat(vegdat.five.Phi(8,:),3,1);
% phi(25:27,:) = repmat(vegdat.five.Phi(9,:),3,1);
% phiw = phi; clear phi
% phi = [phif;phiw]; clear phif phiw
% %extract quarter quadrat statistics as well
% phiqq = zeros(10,1);
% phiqq(1) = vegdat.four.Phi(6,1);
% phiqq(2) = vegdat.four.Phi(7,1);
% phiqq(3) = vegdat.four.Phi(8,1);
% phiqq(4) = vegdat.four.Phi(9,1);
% phiqq(5) = vegdat.four.Phi(10,1);
% phiqq(6) = vegdat.four.Phi(11,1);
% phiqq(7) = vegdat.four.Phi(12,1);
% phiqq(8) = vegdat.four.Phi(13,1);
% phiqq(9) = vegdat.four.Phi(14,1);
% phiqq(10) = vegdat.four.Phi(15,1);

hc = zeros(24,3);
hc(1:3,:) = repmat(vegdat.four.hc(1,:),3,1);
hc(4:6,:) = repmat(vegdat.four.hc(2,:),3,1);
hc(7:9,:) = repmat(vegdat.four.hc(3,:),3,1);
hc(10:12,:) = repmat(vegdat.four.hc(4,:),3,1);
hc(13:15,:) = repmat(vegdat.four.hc(5,:),3,1);
hc(16:18,:) = repmat(vegdat.four.hc(16,:),3,1);
hc(19:21,:) = repmat(vegdat.four.hc(17,:),3,1);
hc(22:24,:) = repmat(vegdat.four.hc(18,:),3,1);
hcf = hc; clear hc
hc = zeros(27,3);
hc(1:3,:) = repmat(vegdat.five.Hc(1,:),3,1);
hc(4:6,:) = repmat(vegdat.five.Hc(2,:),3,1);
hc(7:9,:) = repmat(vegdat.five.Hc(3,:),3,1);
hc(10:12,:) = repmat(vegdat.five.Hc(4,:),3,1);
hc(13:15,:) = repmat(vegdat.five.Hc(5,:),3,1);
hc(16:18,:) = repmat(vegdat.five.Hc(6,:),3,1);
hc(19:21,:) = repmat(vegdat.five.Hc(7,:),3,1);
hc(22:24,:) = repmat(vegdat.five.Hc(8,:),3,1);
hc(25:27,:) = repmat(vegdat.five.Hc(9,:),3,1);
hcw = hc; clear hc
hc = [hcf;hcw]; clear hcf hcw

Z = zeros(24,3);
Z(1:3,:) = repmat(vegdat.four.Z(1,:),3,1);
Z(4:6,:) = repmat(vegdat.four.Z(2,:),3,1);
Z(7:9,:) = repmat(vegdat.four.Z(3,:),3,1);
Z(10:12,:) = repmat(vegdat.four.Z(4,:),3,1);
Z(13:15,:) = repmat(vegdat.four.Z(5,:),3,1);
Z(16:18,:) = repmat(vegdat.four.Z(16,:),3,1);
Z(19:21,:) = repmat(vegdat.four.Z(17,:),3,1);
Z(22:24,:) = repmat(vegdat.four.Z(18,:),3,1);
Zf = Z; clear Z
Z = zeros(27,3);
Z(1:3,:) = repmat(vegdat.five.Z(1,:),3,1);
Z(4:6,:) = repmat(vegdat.five.Z(2,:),3,1);
Z(7:9,:) = repmat(vegdat.five.Z(3,:),3,1);
Z(10:12,:) = repmat(vegdat.five.Z(4,:),3,1);
Z(13:15,:) = repmat(vegdat.five.Z(5,:),3,1);
Z(16:18,:) = repmat(vegdat.five.Z(6,:),3,1);
Z(19:21,:) = repmat(vegdat.five.Z(7,:),3,1);
Z(22:24,:) = repmat(vegdat.five.Z(8,:),3,1);
Z(25:27,:) = repmat(vegdat.five.Z(9,:),3,1);
Zw = Z; clear Z
Z = [Zf;Zw]./1000;clear Zf Zw
zhc = Z./hc; zhc(isinf(zhc)) = 0;

%extract quarter quadrat distances as well
xqq = zeros(10,1);
xqq(1) = vegdat.four.X(6,1);
xqq(2) = vegdat.four.X(7,1);
xqq(3) = vegdat.four.X(8,1);
xqq(4) = vegdat.four.X(9,1);
xqq(5) = vegdat.four.X(10,1);
xqq(6) = vegdat.four.X(11,1);
xqq(7) = vegdat.four.X(12,1);
xqq(8) = vegdat.four.X(13,1);
xqq(9) = vegdat.four.X(14,1);
xqq(10) = vegdat.four.X(15,1);

X = zeros(24,3);
X(1:3,:) = repmat(vegdat.four.X(1,:),3,1);
X(4:6,:) = repmat(vegdat.four.X(2,:),3,1);
X(7:9,:) = repmat(vegdat.four.X(3,:),3,1);
X(10:12,:) = repmat(vegdat.four.X(4,:),3,1);
X(13:15,:) = repmat(vegdat.four.X(5,:),3,1);
X(16:18,:) = repmat(vegdat.four.X(16,:),3,1);
X(19:21,:) = repmat(vegdat.four.X(17,:),3,1);
X(22:24,:) = repmat(vegdat.four.X(18,:),3,1);
Xf = X; clear X
X = zeros(27,3);
X(1:3,:) = repmat(vegdat.five.X(1,:),3,1);
X(4:6,:) = repmat(vegdat.five.X(2,:),3,1);
X(7:9,:) = repmat(vegdat.five.X(3,:),3,1);
X(10:12,:) = repmat(vegdat.five.X(4,:),3,1);
X(13:15,:) = repmat(vegdat.five.X(5,:),3,1);
X(16:18,:) = repmat(vegdat.five.X(6,:),3,1);
X(19:21,:) = repmat(vegdat.five.X(7,:),3,1);
X(22:24,:) = repmat(vegdat.five.X(8,:),3,1);
X(25:27,:) = repmat(vegdat.five.X(9,:),3,1);
Xw = X; clear X
X = [Xf;Xw];clear Xf Xw

eps = [veldat.four.E;veldat.five.E];
avg = [veldat.four.avg;veldat.five.avg];
depth = [veldat.four.depth;veldat.five.depth];
Hs = [veldat.four.wrms;veldat.five.wrms].*sqrt(2);
stdev = [veldat.four.Estd;veldat.five.Estd].*2;

%%%Update: 09/06/2016 use the small-scale phi estimates instead of the 1m^2
%%%estimates to plot figures (can we resolve a pattern with a different
%%%scale of veg density?)
dir1 = 'd:\Projects\Mekong_F2014\DataAnalysis\Paper2\VegStats\';
load([dir1 'VegDat14.mat'])
Phi = NaN(18,3);
Phi(1,:) = vegdat.phi(1);Phi(2,:) = vegdat.phi(2);Phi(3,:) = vegdat.phi(3);Phi(4,:) = 0;
Phi(5,:) = vegdat.phi(4);Phi(6,1) = vegdat.phi(6);Phi(7,1) = vegdat.phi(7);Phi(8,1) = vegdat.phi(8);
Phi(9,1) = vegdat.phi(9);Phi(10,1) = vegdat.phi(10);Phi(11,1) = vegdat.phi(11);Phi(12,1) = vegdat.phi(12);
Phi(13,1) = vegdat.phi(13);Phi(14,1) = vegdat.phi(14);Phi(15,1) = vegdat.phi(5);Phi(16,1) = 0;
Phi(16,2) = vegdat.phi(15);Phi(16,3) = vegdat.phi(16);Phi(17,1) = 0;Phi(17,2) = vegdat.phi(15);
Phi(17,3) = vegdat.phi(16);
Phi(18,:) = vegdat.phi(17);

%extract quarter quadrat statistics as well
phiqq = zeros(10,1);
phiqq(1) = Phi(6,1);
phiqq(2) = Phi(7,1);
phiqq(3) = Phi(8,1);
phiqq(4) = Phi(9,1);
phiqq(5) = Phi(10,1);
phiqq(6) = Phi(11,1);
phiqq(7) = Phi(12,1);
phiqq(8) = Phi(13,1);
phiqq(9) = Phi(14,1);
phiqq(10) = Phi(15,1);

%will need to extract the vegetation data and put it in the same format as
%the velocity data, which is 3 vectrinos as columns, 3*no. of deployments in
%rows (rows 1:3 correspond to 'low' 'mid' and 'high' tide).

phi = zeros(24,3);
phi(1:3,:) = repmat(Phi(1,:),3,1);
phi(4:6,:) = repmat(Phi(2,:),3,1);
phi(7:9,:) = repmat(Phi(3,:),3,1);
phi(10:12,:) = repmat(Phi(4,:),3,1);
phi(13:15,:) = repmat(Phi(5,:),3,1);
phi(16:18,:) = repmat(Phi(16,:),3,1);
phi(19:21,:) = repmat(Phi(17,:),3,1);
phi(22:24,:) = repmat(Phi(18,:),3,1);
phif = phi; clear phi
dir2 = 'd:\Projects\Mekong_W2015\DataAnalysis\Paper2\VegStats\';
load([dir2 'VegDat15.mat'])
Phi = NaN(9,3);
Phi(1,1) = 0;Phi(1,2) = vegdat.phi(1);Phi(1,3) = vegdat.phi(2);
Phi(2,1) = 0;Phi(2,2) = vegdat.phi(1);Phi(2,3) = vegdat.phi(2);
Phi(3,:) = vegdat.phi(3);Phi(4,:) = vegdat.phi(4);Phi(5,:) = vegdat.phi(5);
Phi(6,:) = vegdat.phi(5);
Phi(7,1) = 0;Phi(7,2) = vegdat.phi(6);Phi(7,3) = vegdat.phi(7);
Phi(8,1) = 0;Phi(8,2) = vegdat.phi(6);Phi(8,3) = vegdat.phi(7);
Phi(9,:) = vegdat.phi(8);

phi = zeros(27,3);
phi(1:3,:) = repmat(Phi(1,:),3,1);
phi(4:6,:) = repmat(Phi(2,:),3,1);
phi(7:9,:) = repmat(Phi(3,:),3,1);
phi(10:12,:) = repmat(Phi(4,:),3,1);
phi(13:15,:) = repmat(Phi(5,:),3,1);
phi(16:18,:) = repmat(Phi(6,:),3,1);
phi(19:21,:) = repmat(Phi(7,:),3,1);
phi(22:24,:) = repmat(Phi(8,:),3,1);
phi(25:27,:) = repmat(Phi(9,:),3,1);
phiw = phi; clear phi
phi = [phif;phiw]; clear phif phiw
