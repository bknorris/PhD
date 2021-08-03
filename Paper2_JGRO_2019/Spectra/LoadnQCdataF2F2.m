clear
savedatdir = 'C:\users\bkn5\Projects\Mekong_W2015\DataAnalysis\Paper1\';
tic
disp(['Data loading started at: ' datestr(now)])
%% Load the VPROs
%VPROs from F2F2, Day 2
%Load vecpro files, crop & QC, then save all into a structure

disp('Loading Vectrino Profiler VP1, 06/03/2015')
vpdir = 'C:\Users\bkn5\Projects\Mekong_W2015\Data\Vectrino\6March2015\Vectrinos\';
load([vpdir 'VP1_060315.mat'])
t1 = 7.360295678009259e+05;F2F2.times.t1 = t1;
e1 = 7.360297095254629e+05;F2F2.times.e1 = e1;
VPRO.Data = cropvp(VPRO.Data,t1,e1); %crop data to specified times
[VPRO.Data,VPRO.Config] = VPQC(VPRO.Data,VPRO.Config,2,0); %run QC on data
VPRO.Data = fixbbvpvels(VPRO.Data,2); %VP1-3 on 07/03 need below-bed velocities removed.
%Save important fields into Data structure
F2F2.vpro1.time = VPRO.Data.Time;
F2F2.vpro1.beam1 = VPRO.Data.Profiles_VelBeam1;
F2F2.vpro1.beam2 = VPRO.Data.Profiles_VelBeam2;
F2F2.vpro1.beam3 = VPRO.Data.Profiles_VelBeam3;
F2F2.vpro1.beam4 = VPRO.Data.Profiles_VelBeam4;
[VPRO.Data,VPRO.Config] = VPro_coordinateTransform(VPRO.Data,VPRO.Config,'bx');
F2F2.vpro1.x = VPRO.Data.Profiles_VelX;
F2F2.vpro1.y = VPRO.Data.Profiles_VelY;
F2F2.vpro1.z = (VPRO.Data.Profiles_VelZ1+VPRO.Data.Profiles_VelZ2)./2;
F2F2.vpro1.sr = VPRO.Config.sampleRate;
F2F2.vpro1.rb = VPRO.Data.Profiles_Range;
clear VPRO

disp('Loading Vectrino Profiler VP2, 06/03/2015')
load([vpdir 'VP2_060315.mat'])
VPRO.Data = cropvp(VPRO.Data,t1,e1); %crop data to specified times
[VPRO.Data,VPRO.Config] = VPQC(VPRO.Data,VPRO.Config,2,0); %run QC on data
VPRO.Data = fixbbvpvels(VPRO.Data,2); %VP1-3 on 07/03 need below-bed velocities removed.
%Save important fields into Data structure
F2F2.vpro2.time = VPRO.Data.Time;
F2F2.vpro2.beam1 = VPRO.Data.Profiles_VelBeam1;
F2F2.vpro2.beam2 = VPRO.Data.Profiles_VelBeam2;
F2F2.vpro2.beam3 = VPRO.Data.Profiles_VelBeam3;
F2F2.vpro2.beam4 = VPRO.Data.Profiles_VelBeam4;
[VPRO.Data,VPRO.Config] = VPro_coordinateTransform(VPRO.Data,VPRO.Config,'bx');
F2F2.vpro2.x = VPRO.Data.Profiles_VelX;
F2F2.vpro2.y = VPRO.Data.Profiles_VelY;
F2F2.vpro2.z = (VPRO.Data.Profiles_VelZ1+VPRO.Data.Profiles_VelZ2)./2;
F2F2.vpro2.sr = VPRO.Config.sampleRate;
F2F2.vpro2.rb = VPRO.Data.Profiles_Range;
clear VPRO

disp('Loading Vectrino Profiler VP3, 06/03/2015')
load([vpdir 'VP3_060315.mat'])
VPRO.Data = cropvp(VPRO.Data,t1,e1); %crop data to specified times
[VPRO.Data,VPRO.Config] = VPQC(VPRO.Data,VPRO.Config,2,0); %run QC on data
VPRO.Data = fixbbvpvels(VPRO.Data,2); %VP1-3 on 07/03 need below-bed velocities removed.
%Save important fields into Data structure
F2F2.vpro3.time = VPRO.Data.Time;
F2F2.vpro3.beam1 = VPRO.Data.Profiles_VelBeam1;
F2F2.vpro3.beam2 = VPRO.Data.Profiles_VelBeam2;
F2F2.vpro3.beam3 = VPRO.Data.Profiles_VelBeam3;
F2F2.vpro3.beam4 = VPRO.Data.Profiles_VelBeam4;
[VPRO.Data,VPRO.Config] = VPro_coordinateTransform(VPRO.Data,VPRO.Config,'bx');
F2F2.vpro3.x = VPRO.Data.Profiles_VelX;
F2F2.vpro3.y = VPRO.Data.Profiles_VelY;
F2F2.vpro3.z = (VPRO.Data.Profiles_VelZ1+VPRO.Data.Profiles_VelZ2)./2;
F2F2.vpro3.sr = VPRO.Config.sampleRate;
F2F2.vpro3.rb = VPRO.Data.Profiles_Range;
clear VPRO

disp('Saving F2F2 06/03/2015 velocity file')
save([savedatdir 'F2F2day2Vels'],'F2F2','-v7.3')

%clear F2F2 to save RAM, save times first
times = F2F2.times;
% clearvars -except times gmt2ict savedatdir
% F2F2.times = times;

%clear F2F2 to save RAM
clear

disp(['Data QCd and saved to file in: ' num2str(toc/60) ' minutes'])
