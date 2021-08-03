

% load('D:\Projects\Mekong_F2014\Data\Vectrino\Control\VP1_270914.mat')
% [vp1,~] = VPro_coordinateTransform(VPRO.Data,VPRO.Config,'bx');
% load('D:\Projects\Mekong_F2014\Data\Vectrino\Control\VP3_270914.mat')
% [vp3,~] = VPro_coordinateTransform(VPRO.Data,VPRO.Config,'bx');
% 
% start = datenum(2014,09,27,14,23,57);
% stop = datenum(2014,09,27,14,24,07);
% time1 = vp1.Profiles_HostTimeMatlab-datenum(0,0,0,6,0,0);
% id = find(time1 >= start & time1 <= stop);
% zlw = mean((vp1.Profiles_VelZ1(id,13:17)+vp1.Profiles_VelZ2(id,13:17)./2),2);
% time = time1(id);
% t0 = start;
% 
% 
% time1 = vp3.Profiles_HostTimeMatlab-datenum(0,0,0,6,0,0);
% id = find(time1 >= start & time1 <= stop);
% zup = mean((vp3.Profiles_VelZ1(id,13:17)+vp3.Profiles_VelZ2(id,13:17)./2),2);
% z = [zlw zup];
% time2 = linspace(0,10,length(time1(id)));
% cntrl.datetime = time;
% cntrl.seconds = time2';
% cntrl.z = z;
% cntrl.t0 = t0;
% save('d:\Projects\Mekong_W2015\DataAnalysis\TOS\ControlZVels','cntrl','-v7.3')


load('d:\Projects\Mekong_W2015\Data\Vectrino\14March2015\Vectrinos\VectrinoData.073.15.VP1.00002.mat')
[vp1,~] = VPro_coordinateTransform(Data,Config,'bx');

start = datenum(2015,03,14,05,19,50);
stop = datenum(2015,03,14,05,20,00);
time1 = vp1.Profiles_HostTimeMatlab+datenum(0,0,0,7,0,0);
id = find(time1 >= start & time1 <= stop);
zlw = mean((vp1.Profiles_VelZ1(id,13:17)+vp1.Profiles_VelZ2(id,13:17)./2),2);
time = time1(id);
t0 = start;

load('d:\Projects\Mekong_W2015\Data\Vectrino\14March2015\Vectrinos\VectrinoData.073.15.VP3.00002.mat')
[vp3,~] = VPro_coordinateTransform(Data,Config,'bx');
clear VPRO

time1 = vp3.Profiles_HostTimeMatlab+datenum(0,0,0,7,0,0);
id = find(time1 >= start & time1 <= stop);
zup = mean((vp3.Profiles_VelZ1(id,13:17)+vp3.Profiles_VelZ2(id,13:17)./2),2);
z = [zlw zup];
time2 = linspace(0,10,length(time1(id)));
dps.datetime = time;
dps.seconds = time2';
dps.z = z;
dps.t0 = t0;
save('d:\Projects\Mekong_W2015\DataAnalysis\TOS\DPS2ZVels','dps','-v7.3')