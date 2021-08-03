%try calculating the vorticity. Can only use horizontally-synoptic
%measurments. 
load('D:\Mekong_W2015\DataAnalysis\Paper2\AveragedVelocities\HTA_2Vave.mat')
w1 = (Avgs.vpro2.z1+Avgs.vpro2.z2)./2;
w2 = (Avgs.vpro3.z1+Avgs.vpro3.z2)./2;
v1 = Avgs.vpro2.y;
v2 = Avgs.vpro3.y;
x1 = Avgs.vpro2.x;x2 = Avgs.vpro3.x;
du = nanmean(x2,2);
dv = nanmean(v1-v2,2);
dw = nanmean(w1-w2,2);
% du = nanmean((x1+x2)./2,2);

omgy = du-dw;
omgz = dv-du;