%% Load the data
savedatdir = 'C:\users\bkn5\Projects\Mekong_W2015\DataAnalysis\Spectra\Paper1\QCd\';
tic
disp(['Data loading started at: ' datestr(now)])
%% Load the VPROs
%VPROs from HTA, Day 1
%Load vecpro files, crop & QC, then save all into a structure, we'll call
%it HTA1 for day 1, HTA2 for day 2, and HTA3 for day 3

disp('Loading Vectrino Profiler VP1, 07/03/2015')
vpdir = 'C:\Users\bkn5\Projects\Mekong_W2015\Data\Vectrino\7March2015\Vectrinos\';
load([vpdir 'VP1_070315.mat'])
t1 = HTA.times.t1;e1 = HTA.times.e1;
VPRO.Data = cropvp(VPRO.Data,t1,e1); %crop data to specified times
[VPRO.Data,VPRO.Config] = VPQC(VPRO.Data,VPRO.Config,2,0); %run QC on data
VPRO.Data = fixbbvpvels(VPRO.Data,2); %VP1-3 on 07/03 need below-bed velocities removed.
%Save important fields into Data structure
HTA.vpro1.time = VPRO.Data.Time;
HTA.vpro1.beam1 = VPRO.Data.Profiles_VelBeam1;
HTA.vpro1.beam2 = VPRO.Data.Profiles_VelBeam2;
HTA.vpro1.beam3 = VPRO.Data.Profiles_VelBeam3;
HTA.vpro1.beam4 = VPRO.Data.Profiles_VelBeam4;
[VPRO.Data,VPRO.Config] = VPro_coordinateTransform(VPRO.Data,VPRO.Config,'bx');
HTA.vpro1.x = VPRO.Data.Profiles_VelX;
HTA.vpro1.y = VPRO.Data.Profiles_VelY;
HTA.vpro1.z1 = VPRO.Data.Profiles_VelZ1;
HTA.vpro1.z2 = VPRO.Data.Profiles_VelZ2;
HTA.vpro1.sr = VPRO.Config.sampleRate;
HTA.vpro1.rb = VPRO.Data.Profiles_Range;
clear VPRO

disp('Loading Vectrino Profiler VP2, 07/03/2015')
load([vpdir 'VP2_070315.mat'])
VPRO.Data = cropvp(VPRO.Data,t1,e1); %crop data to specified times
[VPRO.Data,VPRO.Config] = VPQC(VPRO.Data,VPRO.Config,2,0); %run QC on data
VPRO.Data = fixbbvpvels(VPRO.Data,2); %VP1-3 on 07/03 need below-bed velocities removed.
%Save important fields into Data structure
HTA.vpro2.time = VPRO.Data.Time;
HTA.vpro2.beam1 = VPRO.Data.Profiles_VelBeam1;
HTA.vpro2.beam2 = VPRO.Data.Profiles_VelBeam2;
HTA.vpro2.beam3 = VPRO.Data.Profiles_VelBeam3;
HTA.vpro2.beam4 = VPRO.Data.Profiles_VelBeam4;
[VPRO.Data,VPRO.Config] = VPro_coordinateTransform(VPRO.Data,VPRO.Config,'bx');
HTA.vpro2.x = VPRO.Data.Profiles_VelX;
HTA.vpro2.y = VPRO.Data.Profiles_VelY;
HTA.vpro2.z1 = VPRO.Data.Profiles_VelZ1;
HTA.vpro2.z2 = VPRO.Data.Profiles_VelZ2;
HTA.vpro2.sr = VPRO.Config.sampleRate;
HTA.vpro2.rb = VPRO.Data.Profiles_Range;
clear VPRO

disp('Loading Vectrino Profiler VP3, 07/03/2015')
load([vpdir 'VP3_070315.mat'])
VPRO.Data = cropvp(VPRO.Data,t1,e1); %crop data to specified times
[VPRO.Data,VPRO.Config] = VPQC(VPRO.Data,VPRO.Config,2,0); %run QC on data
VPRO.Data = fixbbvpvels(VPRO.Data,2); %VP1-3 on 07/03 need below-bed velocities removed.
%Save important fields into Data structure
HTA.vpro3.time = VPRO.Data.Time;
HTA.vpro3.beam1 = VPRO.Data.Profiles_VelBeam1;
HTA.vpro3.beam2 = VPRO.Data.Profiles_VelBeam2;
HTA.vpro3.beam3 = VPRO.Data.Profiles_VelBeam3;
HTA.vpro3.beam4 = VPRO.Data.Profiles_VelBeam4;
[VPRO.Data,VPRO.Config] = VPro_coordinateTransform(VPRO.Data,VPRO.Config,'bx');
HTA.vpro3.x = VPRO.Data.Profiles_VelX;
HTA.vpro3.y = VPRO.Data.Profiles_VelY;
HTA.vpro3.z1 = VPRO.Data.Profiles_VelZ1;
HTA.vpro3.z2 = VPRO.Data.Profiles_VelZ2;
HTA.vpro3.sr = VPRO.Config.sampleRate;
HTA.vpro3.rb = VPRO.Data.Profiles_Range;
clear VPRO

disp('Saving HTA 07/03/2015 velocity file')
save([savedatdir 'HTAday1Vels'],'HTA','-v7.3')

%clear HTA to save RAM, save times first
times = HTA.times;
clearvars -except times gmt2ict savedatdir VTA
HTA.times = times;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp('Loading Vectrino Profiler VP1, 08/03/2015')
vpdir = 'C:\Users\bkn5\Projects\Mekong_W2015\Data\Vectrino\8March2015\Vectrinos\';
load([vpdir 'VP1_080315.mat'])
t1 = HTA.times.t2;e1 = HTA.times.e2;
VPRO.Data = cropvp(VPRO.Data,t1,e1); %crop data to specified times
[VPRO.Data,VPRO.Config] = VPQC(VPRO.Data,VPRO.Config,2,0); %run QC on data
%Save important fields into Data structure
HTA.vpro1.time = VPRO.Data.Time;
HTA.vpro1.beam1 = VPRO.Data.Profiles_VelBeam1;
HTA.vpro1.beam2 = VPRO.Data.Profiles_VelBeam2;
HTA.vpro1.beam3 = VPRO.Data.Profiles_VelBeam3;
HTA.vpro1.beam4 = VPRO.Data.Profiles_VelBeam4;
[VPRO.Data,VPRO.Config] = VPro_coordinateTransform(VPRO.Data,VPRO.Config,'bx');
HTA.vpro1.x = VPRO.Data.Profiles_VelX;
HTA.vpro1.y = VPRO.Data.Profiles_VelY;
HTA.vpro1.z1 = VPRO.Data.Profiles_VelZ1;
HTA.vpro1.z2 = VPRO.Data.Profiles_VelZ2;
HTA.vpro1.sr = VPRO.Config.sampleRate;
HTA.vpro1.rb = VPRO.Data.Profiles_Range;
clear VPRO

disp('Loading Vectrino Profiler VP2, 08/03/2015')
load([vpdir 'VP2_080315.mat'])
VPRO.Data = cropvp(VPRO.Data,t1,e1); %crop data to specified times
[VPRO.Data,VPRO.Config] = VPQC(VPRO.Data,VPRO.Config,2,0); %run QC on data
%Save important fields into Data structure
HTA.vpro2.time = VPRO.Data.Time;
HTA.vpro2.beam1 = VPRO.Data.Profiles_VelBeam1;
HTA.vpro2.beam2 = VPRO.Data.Profiles_VelBeam2;
HTA.vpro2.beam3 = VPRO.Data.Profiles_VelBeam3;
HTA.vpro2.beam4 = VPRO.Data.Profiles_VelBeam4;
[VPRO.Data,VPRO.Config] = VPro_coordinateTransform(VPRO.Data,VPRO.Config,'bx');
HTA.vpro2.x = VPRO.Data.Profiles_VelX;
HTA.vpro2.y = VPRO.Data.Profiles_VelY;
HTA.vpro2.z1 = VPRO.Data.Profiles_VelZ1;
HTA.vpro2.z2 = VPRO.Data.Profiles_VelZ2;
HTA.vpro2.sr = VPRO.Config.sampleRate;
HTA.vpro2.rb = VPRO.Data.Profiles_Range;
clear VPRO

disp('Loading Vectrino Profiler VP3, 08/03/2015')
load([vpdir 'VP3_080315.mat'])
VPRO.Data = cropvp(VPRO.Data,t1,e1); %crop data to specified times
[VPRO.Data,VPRO.Config] = VPQC(VPRO.Data,VPRO.Config,2,0); %run QC on data
%Save important fields into Data structure
HTA.vpro3.time = VPRO.Data.Time;
HTA.vpro3.beam1 = VPRO.Data.Profiles_VelBeam1;
HTA.vpro3.beam2 = VPRO.Data.Profiles_VelBeam2;
HTA.vpro3.beam3 = VPRO.Data.Profiles_VelBeam3;
HTA.vpro3.beam4 = VPRO.Data.Profiles_VelBeam4;
[VPRO.Data,VPRO.Config] = VPro_coordinateTransform(VPRO.Data,VPRO.Config,'bx');
HTA.vpro3.x = VPRO.Data.Profiles_VelX;
HTA.vpro3.y = VPRO.Data.Profiles_VelY;
HTA.vpro3.z1 = VPRO.Data.Profiles_VelZ1;
HTA.vpro3.z2 = VPRO.Data.Profiles_VelZ2;
HTA.vpro3.sr = VPRO.Config.sampleRate;
HTA.vpro3.rb = VPRO.Data.Profiles_Range;
clear VPRO

disp('Saving HTA 08/03/2015 velocity file')
save([savedatdir 'HTAday2Vels'],'HTA','-v7.3')

%clear HTA to save RAM, save times first
times = HTA.times;
clearvars -except times gmt2ict savedatdir VTA
HTA.times = times;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp('Loading Vectrino Profiler VP1, 10/03/2015')
vpdir = 'C:\Users\bkn5\Projects\Mekong_W2015\Data\Vectrino\10March2015\Vectrinos\';
load([vpdir 'VP1_100315.mat'])
t1 = HTA.times.t3;e1 = HTA.times.e3;
VPRO.Data = cropvp(VPRO.Data,t1,e1); %crop data to specified times
[VPRO.Data,VPRO.Config] = VPQC(VPRO.Data,VPRO.Config,2,0); %run QC on data
%Save important fields into Data structure
HTA.vpro1.time = VPRO.Data.Time;
HTA.vpro1.beam1 = VPRO.Data.Profiles_VelBeam1;
HTA.vpro1.beam2 = VPRO.Data.Profiles_VelBeam2;
HTA.vpro1.beam3 = VPRO.Data.Profiles_VelBeam3;
HTA.vpro1.beam4 = VPRO.Data.Profiles_VelBeam4;
[VPRO.Data,VPRO.Config] = VPro_coordinateTransform(VPRO.Data,VPRO.Config,'bx');
HTA.vpro1.x = VPRO.Data.Profiles_VelX;
HTA.vpro1.y = VPRO.Data.Profiles_VelY;
HTA.vpro1.z1 = VPRO.Data.Profiles_VelZ1;
HTA.vpro1.z2 = VPRO.Data.Profiles_VelZ2;
HTA.vpro1.sr = VPRO.Config.sampleRate;
HTA.vpro1.rb = VPRO.Data.Profiles_Range;
clear VPRO

disp('Loading Vectrino Profiler VP2, 10/03/2015')
load([vpdir 'VP2_100315.mat'])
VPRO.Data = cropvp(VPRO.Data,t1,e1); %crop data to specified times
[VPRO.Data,VPRO.Config] = VPQC(VPRO.Data,VPRO.Config,2,0); %run QC on data
%Save important fields into Data structure
HTA.vpro2.time = VPRO.Data.Time;
HTA.vpro2.beam1 = VPRO.Data.Profiles_VelBeam1;
HTA.vpro2.beam2 = VPRO.Data.Profiles_VelBeam2;
HTA.vpro2.beam3 = VPRO.Data.Profiles_VelBeam3;
HTA.vpro2.beam4 = VPRO.Data.Profiles_VelBeam4;
[VPRO.Data,VPRO.Config] = VPro_coordinateTransform(VPRO.Data,VPRO.Config,'bx');
HTA.vpro2.x = VPRO.Data.Profiles_VelX;
HTA.vpro2.y = VPRO.Data.Profiles_VelY;
HTA.vpro2.z1 = VPRO.Data.Profiles_VelZ1;
HTA.vpro2.z2 = VPRO.Data.Profiles_VelZ2;
HTA.vpro2.sr = VPRO.Config.sampleRate;
HTA.vpro2.rb = VPRO.Data.Profiles_Range;
clear VPRO

disp('Loading Vectrino Profiler VP3, 10/03/2015')
load([vpdir 'VP3_100315.mat'])
VPRO.Data = cropvp(VPRO.Data,t1,e1); %crop data to specified times
[VPRO.Data,VPRO.Config] = VPQC(VPRO.Data,VPRO.Config,2,0); %run QC on data
%Save important fields into Data structure
HTA.vpro3.time = VPRO.Data.Time;
HTA.vpro3.beam1 = VPRO.Data.Profiles_VelBeam1;
HTA.vpro3.beam2 = VPRO.Data.Profiles_VelBeam2;
HTA.vpro3.beam3 = VPRO.Data.Profiles_VelBeam3;
HTA.vpro3.beam4 = VPRO.Data.Profiles_VelBeam4;
[VPRO.Data,VPRO.Config] = VPro_coordinateTransform(VPRO.Data,VPRO.Config,'bx');
HTA.vpro3.x = VPRO.Data.Profiles_VelX;
HTA.vpro3.y = VPRO.Data.Profiles_VelY;
HTA.vpro3.z1 = VPRO.Data.Profiles_VelZ1;
HTA.vpro3.z2 = VPRO.Data.Profiles_VelZ2;
HTA.vpro3.sr = VPRO.Config.sampleRate;
HTA.vpro3.rb = VPRO.Data.Profiles_Range;
clear VPRO

disp('Saving HTA 10/03/2015 velocity file')
save([savedatdir 'HTAday3Vels'],'HTA','-v7.3')

%clear HTA to save RAM
times = HTA.times;
clearvars -except times gmt2ict savedatdir VTA
HTA.times = times;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
datdir2 = 'C:\users\bkn5\Projects\Mekong_W2015\DataAnalysis\Spectra\Paper1\';
ifvta = exist([datdir2 'VTAvelocities' '.mat'],'file');
if ifvta
    %this step takes ages. If it's been run all ready, skip it
elseif ~ifvta
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %VP1 from VTA
    vpdir = 'C:\Users\bkn5\Projects\Mekong_W2015\Data\Vectrino\14March2015\Vectrinos\';
    disp('Loading Vectrino Profiler VP1, 14/03/2015')
    load([vpdir 'VP1_140315.mat'])
    t2 = VTA.times.t2;e2 = VTA.times.e2;
    
    VPRO.Data = cropvp(VPRO.Data,t2,e2); %crop data to specified times
    [VPRO.Data,VPRO.Config] = VPQC(VPRO.Data,VPRO.Config,2,0); %run QC on data
    VPRO.Data = fixbbvpvels(VPRO.Data,2);
    %Save important fields into Data structure  
    VTA.vpro1.time = VPRO.Data.Time;
    VTA.vpro1.beam1 = VPRO.Data.Profiles_VelBeam1;
    VTA.vpro1.beam2 = VPRO.Data.Profiles_VelBeam2;
    VTA.vpro1.beam3 = VPRO.Data.Profiles_VelBeam3;
    VTA.vpro1.beam4 = VPRO.Data.Profiles_VelBeam4;
    [VPRO.Data,VPRO.Config] = VPro_coordinateTransform(VPRO.Data,VPRO.Config,'bx');
    VTA.vpro1.x = VPRO.Data.Profiles_VelX;
    VTA.vpro1.y = VPRO.Data.Profiles_VelY;
    VTA.vpro1.z1 = VPRO.Data.Profiles_VelZ1;
    VTA.vpro1.z2 = VPRO.Data.Profiles_VelZ2;
    VTA.vpro1.sr = VPRO.Config.sampleRate;
    VTA.vpro1.rb = VPRO.Data.Profiles_Range;
    clear VPRO
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %VP2 from VTA
    disp('Loading Vectrino Profiler VP2, 14/03/2015')
    load([vpdir 'VP2_140315.mat'])
    VPRO.Data = cropvp(VPRO.Data,t2,e2); %crop data to specified times
    [VPRO.Data,VPRO.Config] = VPQC(VPRO.Data,VPRO.Config,2,0); %run QC on data
    %Save important fields into Data structure  
    VTA.vpro2.time = VPRO.Data.Time;
    VTA.vpro2.beam1 = VPRO.Data.Profiles_VelBeam1;
    VTA.vpro2.beam2 = VPRO.Data.Profiles_VelBeam2;
    VTA.vpro2.beam3 = VPRO.Data.Profiles_VelBeam3;
    VTA.vpro2.beam4 = VPRO.Data.Profiles_VelBeam4;
    [VPRO.Data,VPRO.Config] = VPro_coordinateTransform(VPRO.Data,VPRO.Config,'bx');
    VTA.vpro2.x = VPRO.Data.Profiles_VelX;
    VTA.vpro2.y = VPRO.Data.Profiles_VelY;
    VTA.vpro2.z1 = VPRO.Data.Profiles_VelZ1;
    VTA.vpro2.z2 = VPRO.Data.Profiles_VelZ2;    
    VTA.vpro2.sr = VPRO.Config.sampleRate;
    VTA.vpro2.rb = VPRO.Data.Profiles_Range;
    clear VPRO   
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %VP3 from VTA
    disp('Loading Vectrino Profiler VP3, 14/03/2015')
    load([vpdir 'VP3_140315.mat'])
    VPRO.Data = cropvp(VPRO.Data,t2,e2); %crop data to specified times
    [VPRO.Data,VPRO.Config] = VPQC(VPRO.Data,VPRO.Config,2,0); %run QC on data
    %Save important fields into Data structure  
    VTA.vpro3.time = VPRO.Data.Time;
    VTA.vpro3.beam1 = VPRO.Data.Profiles_VelBeam1;
    VTA.vpro3.beam2 = VPRO.Data.Profiles_VelBeam2;
    VTA.vpro3.beam3 = VPRO.Data.Profiles_VelBeam3;
    VTA.vpro3.beam4 = VPRO.Data.Profiles_VelBeam4;
    [VPRO.Data,VPRO.Config] = VPro_coordinateTransform(VPRO.Data,VPRO.Config,'bx');
    VTA.vpro3.x = VPRO.Data.Profiles_VelX;
    VTA.vpro3.y = VPRO.Data.Profiles_VelY;
    VTA.vpro3.z1 = VPRO.Data.Profiles_VelZ1;
    VTA.vpro3.z2 = VPRO.Data.Profiles_VelZ2;    
    VTA.vpro3.sr = VPRO.Config.sampleRate;
    VTA.vpro3.rb = VPRO.Data.Profiles_Range;
    clear VPRO
    
    disp('Saving VTA 14/03/2015 velocity file')
    save([datdir2 'VTAvelocities'],'VTA','-v7.3')

    %clear VTA to save RAM
    clearvars -except savedatdir HTA
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end

%% Load the Vectors
%V5109 from HTA, day 1
%I'm going to save the Vectors for the HTA (days 1-3) all in the same file
%to make it easier to load in the future
disp('Loading Vector V5109, 07/03/2015')
vdir = 'C:\Users\bkn5\Projects\Mekong_W2015\Data\Vector\FSS\';
load([vdir 'V5109_070315.mat'])
t1 = HTA.times.t1;e1 = HTA.times.e1;
ind = (find(ADV.datetime >= t1 & ADV.datetime <= e1));
time = ADV.datetime(ind);u = ADV.U(ind);v = ADV.V(ind);w = ADV.W(ind);p = ADV.Pres(ind);
nlin = 1E2;maxg = 1E4;
if cmgidgaps(u)
    u = cmgbridge(u,nlin,maxg,maxg);
end
if cmgidgaps(v)
    v = cmgbridge(v,nlin,maxg,maxg);
end
if cmgidgaps(w)
    w = cmgbridge(w,nlin,maxg,maxg);
end
if cmgidgaps(p)
    p = cmgbridge(w,nlin,maxg,maxg);
end
HTA.day1.V5.time = time;HTA.day1.V5.u = u;HTA.day1.V5.v = v;HTA.day1.V5.w = w;
HTA.day1.V5.p = p;HTA.day1.V5.sr = ADV.Metadata.instmeta.samprate;
clear ADV

disp('Loading Vector V5109, 08/03/2015')
load([vdir 'V5109_070315.mat'])
t1 = HTA.times.t2;e1 = HTA.times.e2;
ind = (find(ADV.datetime >= t1 & ADV.datetime <= e1));
time = ADV.datetime(ind);u = ADV.U(ind);v = ADV.V(ind);w = ADV.W(ind);p = ADV.Pres(ind);
nlin = 1E2;maxg = 1E4;
if cmgidgaps(u)
    u = cmgbridge(u,nlin,maxg,maxg);
end
if cmgidgaps(v)
    v = cmgbridge(v,nlin,maxg,maxg);
end
if cmgidgaps(w)
    w = cmgbridge(w,nlin,maxg,maxg);
end
if cmgidgaps(p)
    p = cmgbridge(w,nlin,maxg,maxg);
end
HTA.day2.V5.time = time;HTA.day2.V5.u = u;HTA.day2.V5.v = v;HTA.day2.V5.w = w;
HTA.day2.V5.p = p;HTA.day2.V5.sr = ADV.Metadata.instmeta.samprate;
clear ADV

disp('Loading Vector V5109, 10/03/2015')
load([vdir 'V5109_100315.mat'])
t1 = HTA.times.t3;e1 = HTA.times.e3;
ind = (find(ADV.datetime >= t1 & ADV.datetime <= e1));
time = ADV.datetime(ind);u = ADV.U(ind);v = ADV.V(ind);w = ADV.W(ind);p = ADV.Pres(ind);
nlin = 1E2;maxg = 1E4;
if cmgidgaps(u)
    u = cmgbridge(u,nlin,maxg,maxg);
end
if cmgidgaps(v)
    v = cmgbridge(v,nlin,maxg,maxg);
end
if cmgidgaps(w)
    w = cmgbridge(w,nlin,maxg,maxg);
end
if cmgidgaps(p)
    p = cmgbridge(w,nlin,maxg,maxg);
end
HTA.day3.V5.time = time;HTA.day3.V5.u = u;HTA.day3.V5.v = v;HTA.day3.V5.w = w;
HTA.day3.V5.p = p;HTA.day3.V5.sr = ADV.Metadata.instmeta.samprate;
clear ADV

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%VC01
disp('Loading Vector VC01, 07/03/2015')
load([vdir 'VC01_070315.mat'])
t1 = HTA.times.t1;e1 = HTA.times.e1;
ind = (find(ADV.datetime >= t1 & ADV.datetime <= e1));
time = ADV.datetime(ind);u = ADV.U(ind);v = ADV.V(ind);w = ADV.W(ind);p = ADV.Pres(ind);
if cmgidgaps(u)
    u = cmgbridge(u,nlin,maxg,maxg);
end
if cmgidgaps(v)
    v = cmgbridge(v,nlin,maxg,maxg);
end
if cmgidgaps(w)
    w = cmgbridge(w,nlin,maxg,maxg);
end
if cmgidgaps(p)
    p = cmgbridge(w,nlin,maxg,maxg);
end
HTA.day1.VC1.time = time;HTA.day1.VC1.u = u;HTA.day1.VC1.v = v;HTA.day1.VC1.w = w;
HTA.day1.VC1.p = p;HTA.day1.VC1.sr = ADV.Metadata.instmeta.samprate;
clear ADV

%VC01
disp('Loading Vector VC01, 08/03/2015')
load([vdir 'VC01_070315.mat'])
t1 = HTA.times.t2;e1 = HTA.times.e2;
ind = (find(ADV.datetime >= t1 & ADV.datetime <= e1));
time = ADV.datetime(ind);u = ADV.U(ind);v = ADV.V(ind);w = ADV.W(ind);p = ADV.Pres(ind);
if cmgidgaps(u)
    u = cmgbridge(u,nlin,maxg,maxg);
end
if cmgidgaps(v)
    v = cmgbridge(v,nlin,maxg,maxg);
end
if cmgidgaps(w)
    w = cmgbridge(w,nlin,maxg,maxg);
end
if cmgidgaps(p)
    p = cmgbridge(w,nlin,maxg,maxg);
end
HTA.day2.VC1.time = time;HTA.day2.VC1.u = u;HTA.day2.VC1.v = v;HTA.day2.VC1.w = w;
HTA.day2.VC1.p = p;HTA.day2.VC1.sr = ADV.Metadata.instmeta.samprate;
clear ADV

disp('Saving HTA Vector velocity file')
save([savedatdir 'HTAVectorVels'],'HTA','-v7.3')

%clear HTA to save RAM
times = HTA.times;
clearvars -except savedatdir times nlin maxg
HTA.times = times;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Load the Aquadopps
%5116, day 1. Save Aquadopp files as one structure for the ease of loading
disp('Loading Aquadopp AD5116, 07/03/2015')
adir = 'C:\Users\bkn5\Projects\Mekong_W2015\Data\Aquadopp\FSS\';
load([adir 'AD5116_9March2015_f_pad.mat'])
t1 = HTA.times.t1;e1 = HTA.times.e1;
ind = (find(aqdp.datenum >= t1 & aqdp.datenum <= e1));
time = aqdp.datenum(ind,:);u = aqdp.u(ind,:);v = aqdp.v(ind,:);w = aqdp.w(ind,:);p = aqdp.pressure(ind);
if cmgidgaps(u)
    u = cmgbridge(u,nlin,maxg,maxg);
end
if cmgidgaps(v)
    v = cmgbridge(v,nlin,maxg,maxg);
end
if cmgidgaps(w)
    w = cmgbridge(w,nlin,maxg,maxg);
end
if cmgidgaps(p)
    p = cmgbridge(w,nlin,maxg,maxg);
end
HTA.day1.ad6.time = time;HTA.day1.ad6.u = u;HTA.day1.ad6.v = v;HTA.day1.ad6.w = w;HTA.day1.ad6.sr = 8;
HTA.day1.ad6.p = p;hab = aqdp.metadata.HAB/1000;HTA.day1.ad6.rbdepth = hab-aqdp.rangebins;
clear aqdp

disp('Loading Aquadopp AD5116, 08/03/2015')
load([adir 'AD5116_9March2015_f_pad.mat'])
t1 = HTA.times.t2;e1 = HTA.times.e2;
ind = (find(aqdp.datenum >= t1 & aqdp.datenum <= e1));
time = aqdp.datenum(ind,:);u = aqdp.u(ind,:);v = aqdp.v(ind,:);w = aqdp.w(ind,:);p = aqdp.pressure(ind);
if cmgidgaps(u)
    u = cmgbridge(u,nlin,maxg,maxg);
end
if cmgidgaps(v)
    v = cmgbridge(v,nlin,maxg,maxg);
end
if cmgidgaps(w)
    w = cmgbridge(w,nlin,maxg,maxg);
end
if cmgidgaps(p)
    p = cmgbridge(w,nlin,maxg,maxg);
end
HTA.day2.ad6.time = time;HTA.day2.ad6.u = u;HTA.day2.ad6.v = v;HTA.day2.ad6.w = w;HTA.day2.ad6.sr = 8;
HTA.day2.ad6.p = p;hab = aqdp.metadata.HAB/1000;HTA.day2.ad6.rbdepth = hab-aqdp.rangebins;
clear aqdp

disp('Loading Aquadopp AD5116, 10/03/2015')
load([adir 'AD5116_12March2015_f_pad.mat'])
t1 = HTA.times.t3;e1 = HTA.times.e3;
ind = (find(aqdp.datenum >= t1 & aqdp.datenum <= e1));
time = aqdp.datenum(ind,:);u = aqdp.u(ind,:);v = aqdp.v(ind,:);w = aqdp.w(ind,:);p = aqdp.pressure(ind);
if cmgidgaps(u)
    u = cmgbridge(u,nlin,maxg,maxg);
end
if cmgidgaps(v)
    v = cmgbridge(v,nlin,maxg,maxg);
end
if cmgidgaps(w)
    w = cmgbridge(w,nlin,maxg,maxg);
end
if cmgidgaps(p)
    p = cmgbridge(w,nlin,maxg,maxg);
end
HTA.day3.ad6.time = time;HTA.day3.ad6.u = u;HTA.day3.ad6.v = v;HTA.day3.ad6.w = w;HTA.day3.ad6.sr = 8;
HTA.day3.ad6.p = p;hab = aqdp.metadata.HAB/1000;HTA.day3.ad6.rbdepth = hab-aqdp.rangebins;
clear aqdp

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%AD5117, day 1
disp('Loading Aquadopp AD5117, 07/03/2015')
adir = 'C:\Users\bkn5\Projects\Mekong_W2015\Data\Aquadopp\FSS\';
load([adir 'AD5117_9March2015_f_pad.mat'])
t1 = HTA.times.t1;e1 = HTA.times.e1;
ind = (find(aqdp.datenum >= t1 & aqdp.datenum <= e1));
time = aqdp.datenum(ind,:);u = aqdp.u(ind,:);v = aqdp.v(ind,:);w = aqdp.w(ind,:);p = aqdp.pressure(ind);
if cmgidgaps(u)
    u = cmgbridge(u,nlin,maxg,maxg);
end
if cmgidgaps(v)
    v = cmgbridge(v,nlin,maxg,maxg);
end
if cmgidgaps(w)
    w = cmgbridge(w,nlin,maxg,maxg);
end
if cmgidgaps(p)
    p = cmgbridge(w,nlin,maxg,maxg);
end
HTA.day1.ad7.time = time;HTA.day1.ad7.u = u;HTA.day1.ad7.v = v;HTA.day1.ad7.w = w;HTA.day1.ad7.sr = 8;
HTA.day1.ad7.p = p;hab = aqdp.metadata.HAB/1000;HTA.day1.ad7.rbdepth = hab-aqdp.rangebins;
clear aqdp

disp('Loading Aquadopp AD5117, 08/03/2015')
load([adir 'AD5117_9March2015_f_pad.mat'])
t1 = HTA.times.t2;e1 = HTA.times.e2;
ind = (find(aqdp.datenum >= t1 & aqdp.datenum <= e1));
time = aqdp.datenum(ind,:);u = aqdp.u(ind,:);v = aqdp.v(ind,:);w = aqdp.w(ind,:);p = aqdp.pressure(ind);
if cmgidgaps(u)
    u = cmgbridge(u,nlin,maxg,maxg);
end
if cmgidgaps(v)
    v = cmgbridge(v,nlin,maxg,maxg);
end
if cmgidgaps(w)
    w = cmgbridge(w,nlin,maxg,maxg);
end
if cmgidgaps(p)
    p = cmgbridge(w,nlin,maxg,maxg);
end
HTA.day2.ad7.time = time;HTA.day2.ad7.u = u;HTA.day2.ad7.v = v;HTA.day2.ad7.w = w;HTA.day2.ad7.sr = 8;
HTA.day2.ad7.p = p;hab = aqdp.metadata.HAB/1000;HTA.day2.ad7.rbdepth = hab-aqdp.rangebins;
clear aqdp

disp('Loading Aquadopp AD5117, 10/03/2015')
adir = 'C:\Users\bkn5\Projects\Mekong_W2015\Data\Aquadopp\FSS\';
load([adir 'AD5117_12March2015_f_pad.mat'])
t1 = HTA.times.t3;e1 = HTA.times.e3;
ind = (find(aqdp.datenum >= t1 & aqdp.datenum <= e1));
time = aqdp.datenum(ind,:);u = aqdp.u(ind,:);v = aqdp.v(ind,:);w = aqdp.w(ind,:);p = aqdp.pressure(ind);
if cmgidgaps(u)
    u = cmgbridge(u,nlin,maxg,maxg);
end
if cmgidgaps(v)
    v = cmgbridge(v,nlin,maxg,maxg);
end
if cmgidgaps(w)
    w = cmgbridge(w,nlin,maxg,maxg);
end
if cmgidgaps(p)
    p = cmgbridge(w,nlin,maxg,maxg);
end
HTA.day3.ad7.time = time;HTA.day3.ad7.u = u;HTA.day3.ad7.v = v;HTA.day3.ad7.w = w;HTA.day3.ad7.sr = 8;
HTA.day3.ad7.p = p;hab = aqdp.metadata.HAB/1000;HTA.day3.ad7.rbdepth = hab-aqdp.rangebins;
clear aqdp

disp('Saving HTA Aquadopp velocity file')
save([savedatdir 'HTAaqdpVels'],'HTA','-v7.3')

%clear HTA to save RAM
clear

disp(['Data QCd and saved to file in: ' num2str(toc/60) ' minutes'])

