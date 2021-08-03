%Plot script for vertical spectra plots. In this routine, I will (for now)
%include instruments from the Fine Scale Study (2) (now renamed Horizontal
%Turbulence Study, HTA in Paper 1): VP3, V5109 and VC01. I will also
%include AD5116 and AD5117 from this experiment, as well as VP1-VP3 from
%the Dense Pneumatophore Study (2) (now renamed Vertical Turbulence Study,
%VTA). These plots will be gridded on a 3x3 subplot array showing vertical
%flux of water through spectral analysis. Simply put, this script loads the
%data, bridges gaps, computes Welsh's power spectral density function, and
%plots the results.

clear
close all
Data = struct();HTA = struct();VTA = struct();
skipPlots = 0;
skipVTA = 0;
skipSpectra = 0;
%%
%Designate a time cap. This should correspond to the limits where the top
%instrument was in the water. For reference, it will be helpful to have a
%pressure signal plotted with this time limit designated. Load data from
%V5108 (mudflat) and plot the pressure signal of 08/03/2015 with the
%section designated here greyed out.
vdir = 'C:\Users\bkn5\Projects\Mekong_W2015\Data\Vector\SW_Mudflat\';
t1 = datenum(2015,03,08,13,40,00);e1 = datenum(2015,03,08,19,25,00);
%t2 and e2 are the times when the top Vector of the HTA is underwater
t2 = datenum(2015,03,08,15,20,00);e2 = datenum(2015,03,08,17,10,00);
Data.hta.times.t1 = t1;Data.hta.times.e1 = e1;Data.hta.times.t2 = t2;Data.hta.times.e2 = e2;
timesamples = [7.360316435908566e+05;7.360316790364584e+05;7.360317076099537e+05];
Data.hta.times.samples = timesamples;

if ~skipPlots
    load([vdir 'V5108_080315.mat'])
    ind = find(ADV.datetime >= t1 & ADV.datetime <= e1);
    advtime = downsample(ADV.datetime(ind),1000);advpres = downsample(ADV.Pres(ind),1000);
    advdepth = advpres+(ADV.Metadata.pressure_sensor_height/1000);
    sid = zeros(3,1);
    for i = 1:3
        sid(i) = find(advtime == timesamples(i));
    end
    sampledepth = advdepth(sid);
    %timesamples denote the timing of when the spectra will be run. 3 Samples,
    %one rising tide, one at high tide, one falling tide. These were taken from
    %the plot below:
    
    %plot data range figure
    f1 = figure(1);
    set(f1,'PaperOrientation','portrait',...
    'position',[400 100   800 600]);
    set(gcf,'color','w','PaperPositionMode','auto')
    hold on
    c = [0.7 0.7 0.7;0.5 0.5 0.5;0 0 0];
    p = patch([t2 e2 e2 t2],[0 0 2 2],[.95 .95 .95]);set(p,'EdgeColor','none')
    plot(advtime,advdepth,'k','LineWidth',1)
    for i = 1:3
        plot(timesamples(i),sampledepth(i),'o','MarkerSize',10,'MarkerEdgeColor',c(i,:),'MarkerFaceColor',c(i,:),'LineWidth',1.5)
    end
    hold off
    set(gca,'Ylim',[0 1.6],'Xlim',[t1 e1],'box','on','LineWidth',1.5,'FontSize',12)
    datetick('x','HH:MM:SS','keepticks','keeplimits')
    xlabel('\bf\itTime on 08/03/2015')
    ylabel('\bf\itDepth (m)'), hold off
    clear ADV
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Now, do the same for the mudflat Vector from the VTA
vdir = 'C:\Users\bkn5\Projects\Mekong_W2015\Data\Vector\DPS2\';
t1 = datenum(2015,03,14,04,10,00);e1 = datenum(2015,03,14,15,00,00);
t2 = datenum(2015,03,14,07,00,00);e2 = datenum(2015,03,14,12,50,00);
%t2 and e2 are the times when the top Vectrino of the VTA is underwater
Data.vta.times.t1 = t1;Data.vta.times.e1 = e1;Data.vta.times.t2 = t2;Data.vta.times.e2 = e2;
timesamples = [7.360373002025463e+05;7.360373833912037e+05;7.360374886429397e+05];
Data.vta.times.samples = timesamples;

if ~skipPlots
    load([vdir 'V5109_130315.mat'])
    ind = find(ADV.datetime >= t1 & ADV.datetime <= e1);
    advtime = downsample(ADV.datetime(ind),1000);advpres = downsample(ADV.Pres(ind),1000);
    advdepth = advpres+(ADV.Metadata.pressure_sensor_height/1000);
    sid = zeros(3,1);
    for i = 1:3
        sid(i) = find(advtime == timesamples(i));
    end
    sampledepth = advdepth(sid);
    %timesamples denote the timing of when the spectra will be run. 3 Samples,
    %one rising tide, one at high tide, one falling tide. These were taken from
    %the plot below:

    %plot data range figure
    f2 = figure(2);
    set(f2,'PaperOrientation','portrait',...
    'position',[400 100   800 600]);
    set(gcf,'color','w','PaperPositionMode','auto')
    hold on
    c = [0.7 0.7 0.7;0.5 0.5 0.5;0 0 0];
    p = patch([t2 e2 e2 t2],[0 0 2 2],[.95 .95 .95]);set(p,'EdgeColor','none')
    plot(advtime,advdepth,'k','LineWidth',1)
    for i = 1:3
        plot(timesamples(i),sampledepth(i),'o','MarkerSize',10,'MarkerEdgeColor',c(i,:),'MarkerFaceColor',c(i,:),'LineWidth',1.5)
    end
    hold off
    set(gca,'Ylim',[0 1.8],'Xlim',[t1 e1],'box','on','LineWidth',1.5,'FontSize',12)
    datetick('x','HH:MM:SS','keepticks','keeplimits')
    xlabel('\bf\itTime on 14/03/2015')
    ylabel('\bf\itDepth (m)'), hold off
    clear ADV
end

%% Load the data
%Check first to see if the data file has all ready been created. Otherwise,
%skip the loading step and instead load the data file
savedatdir = 'C:\users\bkn5\Projects\Mekong_W2015\DataAnalysis\Spectra\Paper1\';
ifvta = exist([savedatdir 'VTAvelocities' '.mat'],'file');
ifhta = exist([savedatdir 'HTAvelocities' '.mat'],'file');
if ~ifvta || ~ifhta
    tic
    disp(['Data loading started at: ' datestr(now)])
    % Load the VPROs
    %VPRO from HTA
    vpdir = 'C:\Users\bkn5\Projects\Mekong_W2015\Data\Vectrino\8March2015\Vectrinos\';
    load([vpdir 'VP3_080315.mat'])
    gmt2ict = datenum(0,0,0,1,0,0)*7;
    t2 = Data.hta.times.t2;e2 = Data.hta.times.e2;
    VPRO.Data = cropvp(VPRO.Data,t2,e2); %crop data to specified times
    [VPRO.Data,VPRO.Config] = VPQC(VPRO.Data,VPRO.Config,0); %run QC on data
    [VPRO.Data,VPRO.Config] = VPro_coordinateTransform(VPRO.Data,VPRO.Config,'bx');
    
    %Save important fields into Data structure
    Data.hta.vpro.time = VPRO.Data.Time+gmt2ict;
    Data.hta.vpro.x = VPRO.Data.Profiles_VelX;
    Data.hta.vpro.y = VPRO.Data.Profiles_VelY;
    Data.hta.vpro.z = (VPRO.Data.Profiles_VelZ1+VPRO.Data.Profiles_VelZ2)./2;
    Data.hta.vpro.sr = VPRO.Config.sampleRate;
    clear VPRO
    
    if skipVTA || ifvta
        %this step takes ages. If it's been run all ready, skip it
    elseif ~skipVTA && ~ifvta
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %VP1 from VTA
        vpdir = 'C:\Users\bkn5\Projects\Mekong_W2015\Data\Vectrino\14March2015\Vectrinos\';
        load([vpdir 'VP1_140315.mat'])
        t2 = Data.vta.times.t2;e2 = Data.vta.times.e2;
        gmt2ict = datenum(0,0,0,1,0,0)*7;
        VPRO.Data = cropvp(VPRO.Data,t2,e2); %crop data to specified times
        [VPRO.Data,VPRO.Config] = VPQC(VPRO.Data,VPRO.Config,0); %run QC on data
        [VPRO.Data,VPRO.Config] = VPro_coordinateTransform(VPRO.Data,VPRO.Config,'bx');
        
        %Save important fields into Data structure
        Data.vta.vpro1.time = VPRO.Data.Time+gmt2ict;
        Data.vta.vpro1.x = VPRO.Data.Profiles_VelX;
        Data.vta.vpro1.y = VPRO.Data.Profiles_VelY;
        Data.vta.vpro1.z = (VPRO.Data.Profiles_VelZ1+VPRO.Data.Profiles_VelZ2)./2;
        Data.vta.vpro1.sr = VPRO.Config.sampleRate;
        clear VPRO
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %VP2 from VTA
        load([vpdir 'VP2_140315.mat'])
        gmt2ict = datenum(0,0,0,1,0,0)*7;
        VPRO.Data = cropvp(VPRO.Data,t2,e2); %crop data to specified times
        [VPRO.Data,VPRO.Config] = VPQC(VPRO.Data,VPRO.Config,0); %run QC on data
        [VPRO.Data,VPRO.Config] = VPro_coordinateTransform(VPRO.Data,VPRO.Config,'bx');
        
        %Save important fields into Data structure
        Data.vta.vpro2.time = VPRO.Data.Time+gmt2ict;
        Data.vta.vpro2.x = VPRO.Data.Profiles_VelX;
        Data.vta.vpro2.y = VPRO.Data.Profiles_VelY;
        Data.vta.vpro2.z = (VPRO.Data.Profiles_VelZ1+VPRO.Data.Profiles_VelZ2)./2;
        Data.vta.vpro2.sr = VPRO.Config.sampleRate;
        clear VPRO
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %VP3 from VTA
        load([vpdir 'VP3_140315.mat'])
        gmt2ict = datenum(0,0,0,1,0,0)*7;
        VPRO.Data = cropvp(VPRO.Data,t2,e2); %crop data to specified times
        [VPRO.Data,VPRO.Config] = VPQC(VPRO.Data,VPRO.Config,0); %run QC on data
        [VPRO.Data,VPRO.Config] = VPro_coordinateTransform(VPRO.Data,VPRO.Config,'bx');
        
        %Save important fields into Data structure
        Data.vta.vpro3.time = VPRO.Data.Time+gmt2ict;
        Data.vta.vpro3.x = VPRO.Data.Profiles_VelX;
        Data.vta.vpro3.y = VPRO.Data.Profiles_VelY;
        Data.vta.vpro3.z = (VPRO.Data.Profiles_VelZ1+VPRO.Data.Profiles_VelZ2)./2;
        Data.vta.vpro3.sr = VPRO.Config.sampleRate;
        clear VPRO
    end
    % Load the Vectors
    %V5109
    vdir = 'C:\Users\bkn5\Projects\Mekong_W2015\Data\Vector\FSS\';
    load([vdir 'V5109_070315.mat'])
    t2 = Data.hta.times.t2;e2 = Data.hta.times.e2;
    ind = (find(ADV.datetime >= t2 & ADV.datetime <= e2));
    time = ADV.datetime(ind);u = ADV.U(ind);v = ADV.V(ind);w = ADV.W(ind);
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
    Data.hta.midvec.time = time;Data.hta.midvec.u = u;Data.hta.midvec.v = v;Data.hta.midvec.w = w;
    Data.hta.midvec.sr = ADV.Metadata.instmeta.samprate;
    clear ADV
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %VC01
    load([vdir 'VC01_070315.mat'])
    ind = (find(ADV.datetime >= t2 & ADV.datetime <= e2));
    time = ADV.datetime(ind);u = ADV.U(ind);v = ADV.V(ind);w = ADV.W(ind);
    if cmgidgaps(u)
        u = cmgbridge(u,nlin,maxg,maxg);
    end
    if cmgidgaps(v)
        v = cmgbridge(v,nlin,maxg,maxg);
    end
    if cmgidgaps(w)
        w = cmgbridge(w,nlin,maxg,maxg);
    end
    Data.hta.topvec.time = time;Data.hta.topvec.u = u;Data.hta.topvec.v = v;Data.hta.topvec.w = w;
    Data.hta.topvec.sr = ADV.Metadata.instmeta.samprate;
    clear ADV
    
    % Load the Aquadopps
    %5116
    adir = 'C:\Users\bkn5\Projects\Mekong_W2015\Data\Aquadopp\FSS\';
    load([adir 'AD5116_9March2015_f_pad.mat'])
    ind = (find(aqdp.datenum >= t2 & aqdp.datenum <= e2));
    time = aqdp.datenum(ind,:);u = aqdp.u(ind,:);v = aqdp.v(ind,:);w = aqdp.w(ind,:);
    if cmgidgaps(u)
        u = cmgbridge(u,nlin,maxg,maxg);
    end
    if cmgidgaps(v)
        v = cmgbridge(v,nlin,maxg,maxg);
    end
    if cmgidgaps(w)
        w = cmgbridge(w,nlin,maxg,maxg);
    end
    Data.hta.ad6.time = time;Data.hta.ad6.u = u;Data.hta.ad6.v = v;Data.hta.ad6.w = w;Data.hta.ad6.sr = 8;
    hab = aqdp.metadata.HAB/1000;Data.hta.ad6.rbdepth = hab-aqdp.rangebins;
    clear aqdp
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %AD5117
    adir = 'C:\Users\bkn5\Projects\Mekong_W2015\Data\Aquadopp\FSS\';
    load([adir 'AD5117_9March2015_f_pad.mat'])
    ind = (find(aqdp.datenum >= t2 & aqdp.datenum <= e2));
    time = aqdp.datenum(ind,:);u = aqdp.u(ind,:);v = aqdp.v(ind,:);w = aqdp.w(ind,:);
    if cmgidgaps(u)
        u = cmgbridge(u,nlin,maxg,maxg);
    end
    if cmgidgaps(v)
        v = cmgbridge(v,nlin,maxg,maxg);
    end
    if cmgidgaps(w)
        w = cmgbridge(w,nlin,maxg,maxg);
    end
    Data.hta.ad7.time = time;Data.hta.ad7.u = u;Data.hta.ad7.v = v;Data.hta.ad7.w = w;Data.hta.ad7.sr = 8;
    hab = aqdp.metadata.HAB/1000;Data.hta.ad7.rbdepth = hab-aqdp.rangebins;
    clear aqdp
    
    %I'm going to split up the save files since they're HUGE! It's that
    %damn 50Hz Vectrino Data
    VTA = Data.vta;HTA = Data.hta;
    if ~ifvta
        disp('Saving VTA velocity file')
        save([savedatdir 'VTAvelocities'],'VTA','-v7.3') %note the VTA vecpro file is ~10x larger than the HTA file
    end
    if ~ifhta
        disp('Saving HTA velocity file')
        save([savedatdir 'HTAvelocities'],'HTA','-v7.3')
    end
    disp(['Data loaded in: ' num2str(toc/60) ' minutes'])
end
if ifvta
    disp('Loading VTA file')
    load([savedatdir 'VTAvelocities.mat'])
end
if ifhta
    disp('Loading HTA file')
    load([savedatdir 'HTAvelocities.mat'])
end
if ~skipSpectra
    disp('User requests to re-run spectral analysis')
    if isfield(HTA.vpro,'Sz')
        %remove fields so spectral analysis can be rerun
        fn = fieldnames(HTA);
        for i = 2:length(fn)
            HTA.(fn{i}) = rmfield(HTA.(fn{i}),'f');
            if isfield(HTA.(fn{i}),'Sz')
                HTA.(fn{i}) = rmfield(HTA.(fn{i}),'Sz');
            elseif isfield(HTA.(fn{i}),'Sz1')
                fn2 = fieldnames(HTA.(fn{i}));
                for ii = 7:9
                    HTA.(fn{i}) = rmfield(HTA.(fn{i}),fn2{ii});
                end
            end
        end
    end
    if isfield(VTA.vpro1,'Sz')
        fn = fieldnames(VTA);
        for i = 2:length(fn)
            VTA.(fn{i}) = rmfield(VTA.(fn{i}),'Sz');
            VTA.(fn{i}) = rmfield(VTA.(fn{i}),'f');
        end
    end
    %% Perform spectral analysis, compute PSD of vertical (z or w) timeseries
    intv = 10; %minutes
    %Need to select which times to use of the whole t-s
    %VPRO from hta
    fs = HTA.vpro.sr;sr = datenum(0,0,0,0,0,1/fs);
    samples = Data.hta.times.samples;
    sid  = zeros(1,3);
    for i = 1:3
       sdx = find(abs(HTA.vpro.time-samples(i))<=sr); 
       if length(sdx) > 1
           sid(i) = sdx(1);
       else
           sid(i) = sdx;
       end
    end
    avt = fs*intv*60; %#of samples in interval
    minf = 0.005;maxf = 10; %Hz
    idx = [sid(1) sid(1)+avt; sid(2) sid(2)+avt; sid(3) sid(3)+avt];
    for i = 1:length(idx)
        ind = idx(i,1):idx(i,2);
        z = HTA.vpro.z(ind,15); %use bin 15 for vecpros
        if rem(length(z),2) == 1
            n = length(z)-1; %force inputs to be even
            z = z(1:end-1);
        else
            n = length(z);
        end
        z = detrend(z);
        nfft = n/2;
        window = hanning(n);
        noverlap = nfft/2;
        [Sz,f] = pwelch(z,window,noverlap,nfft,fs); %compute PSD
        fr = find(f >= minf & f <= maxf); %set frequency cutoff
        Sz = Sz(fr);f = f(fr);
        HTA.vpro.Sz(:,i) = downsample(Sz,10);HTA.vpro.f(:,i) = downsample(f,10);
    end
    
    %VP1 from VTA
    samples = Data.vta.times.samples;
    sid  = zeros(1,3);
    for i = 1:3
       sdx = find(abs(VTA.vpro1.time-samples(i))<=sr); 
       if length(sdx) > 1
           sid(i) = sdx(1);
       else
           sid(i) = sdx;
       end
    end
    idx = [sid(1) sid(1)+avt; sid(2) sid(2)+avt; sid(3) sid(3)+avt];
    for i = 1:length(idx)
        ind = idx(i,1):idx(i,2);
        z = VTA.vpro1.z(ind,15); %use bin 15 for vecpros
        if rem(length(z),2) == 1
            n = length(z)-1; %force inputs to be even
            z = z(1:end-1);
        else
            n = length(z);
        end
        z = detrend(z);
        nfft = n/2;
        window = hanning(n);
        noverlap = nfft/2;
        [Sz,f] = pwelch(z,window,noverlap,nfft,fs); %compute PSD
        fr = find(f >= minf & f <= maxf); %set frequency cutoff
        Sz = Sz(fr);f = f(fr);
        VTA.vpro1.Sz(:,i) = downsample(Sz,10);VTA.vpro1.f(:,i) = downsample(f,10);
    end
    
    %VP2 from VTA
    sid  = zeros(1,3);
    for i = 1:3
       sdx = find(abs(VTA.vpro2.time-samples(i))<=sr); 
       if length(sdx) > 1
           sid(i) = sdx(1);
       else
           sid(i) = sdx;
       end
    end
    idx = [sid(1) sid(1)+avt; sid(2) sid(2)+avt; sid(3) sid(3)+avt];
    for i = 1:length(idx)
        ind = idx(i,1):idx(i,2);
        z = VTA.vpro2.z(ind,15); %use bin 15 for vecpros
        if rem(length(z),2) == 1
            n = length(z)-1; %force inputs to be even
            z = z(1:end-1);
        else
            n = length(z);
        end
        z = detrend(z);
        nfft = n/2;
        window = hanning(n);
        noverlap = nfft/2;
        [Sz,f] = pwelch(z,window,noverlap,nfft,fs); %compute PSD
        fr = find(f >= minf & f <= maxf); %set frequency cutoff
        Sz = Sz(fr);f = f(fr);
        VTA.vpro2.Sz(:,i) = downsample(Sz,10);VTA.vpro2.f(:,i) = downsample(f,10);
    end
    
    %VP3 from VTA
    sid  = zeros(1,3);
    for i = 1:3
       sdx = find(abs(VTA.vpro3.time-samples(i))<=sr); 
       if length(sdx) > 1
           sid(i) = sdx(1);
       else
           sid(i) = sdx;
       end
    end
    idx = [sid(1) sid(1)+avt; sid(2) sid(2)+avt; sid(3) sid(3)+avt];
    for i = 1:length(idx)
        ind = idx(i,1):idx(i,2);
        z = VTA.vpro3.z(ind,15); %use bin 15 for vecpros
        if rem(length(z),2) == 1
            n = length(z)-1; %force inputs to be even
            z = z(1:end-1);
        else
            n = length(z);
        end
        z = detrend(z);
        nfft = n/2;
        window = hanning(n);
        noverlap = nfft/2;
        [Sz,f] = pwelch(z,window,noverlap,nfft,fs); %compute PSD
        fr = find(f >= minf & f <= maxf); %set frequency cutoff
        Sz = Sz(fr);f = f(fr);
        VTA.vpro3.Sz(:,i) = downsample(Sz,10);VTA.vpro3.f(:,i) = downsample(f,10);
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Vectors from HTA
    %V5109
    fs = HTA.midvec.sr;sr = datenum(0,0,0,0,0,1/fs);
    samples = Data.hta.times.samples;
    sid  = zeros(1,3);
    for i = 1:3
       sdx = find(abs(HTA.midvec.time-samples(i))<=sr); 
       if length(sdx) > 1
           sid(i) = sdx(1);
       else
           sid(i) = sdx;
       end
    end
    avt = fs*intv*60; %#of samples in interval
    minf = 0.005;maxf = 10; %Hz
    idx = [sid(1) sid(1)+avt; sid(2) sid(2)+avt; sid(3) sid(3)+avt];
    for i = 1:length(idx)
        ind = idx(i,1):idx(i,2);
        w = HTA.midvec.w(ind); %use bin 15 for vecpros
        if rem(length(w),2) == 1
            n = length(w)-1; %force inputs to be even
            w = w(1:end-1);
        else
            n = length(w);
        end
        w = detrend(w);
        nfft = n/2;
        window = hanning(n);
        noverlap = nfft/2;
        [Sz,f] = pwelch(w,window,noverlap,nfft,fs); %compute PSD
        fr = find(f >= minf & f <= maxf); %set frequency cutoff
        Sz = Sz(fr);f = f(fr);
        HTA.midvec.Sz(:,i) = downsample(Sz,10);HTA.midvec.f(:,i) = downsample(f,10);
    end
    
    %VC01
    sid  = zeros(1,3);
    for i = 1:3
       sdx = find(abs(HTA.topvec.time-samples(i))<=sr); 
       if length(sdx) > 1
           sid(i) = sdx(1);
       else
           sid(i) = sdx;
       end
    end
    idx = [sid(1) sid(1)+avt; sid(2) sid(2)+avt; sid(3) sid(3)+avt];
    for i = 1:length(idx)
        ind = idx(i,1):idx(i,2);
        w = HTA.topvec.w(ind); %use bin 15 for vecpros
        if rem(length(w),2) == 1
            n = length(w)-1; %force inputs to be even
            w = w(1:end-1);
        else
            n = length(w);
        end
        w = detrend(w);
        nfft = n/2;
        window = hanning(n);
        noverlap = nfft/2;
        [Sz,f] = pwelch(w,window,noverlap,nfft,fs); %compute PSD
        fr = find(f >= minf & f <= maxf); %set frequency cutoff
        Sz = Sz(fr);f = f(fr);
        HTA.topvec.Sz(:,i) = downsample(Sz,10);HTA.topvec.f(:,i) = downsample(f,10);
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Aquadopps from HTA
    %AD5116
    bin = [1 8 16]; %bins to plot spectra
    fs = HTA.ad6.sr;sr = datenum(0,0,0,0,0,1/fs);
    sid  = zeros(1,3);
    for i = 1:3
       sdx = find(abs(HTA.ad6.time-samples(i))<=sr); 
       if length(sdx) > 1
           sid(i) = sdx(1);
       else
           sid(i) = sdx;
       end
    end
    avt = fs*intv*60; %#of samples in interval
    minf = 0.005;maxf = 10; %Hz
    idx = [sid(1) sid(1)+avt; sid(2) sid(2)+avt; sid(3) sid(3)+avt];
    for i = 1:length(idx)
        ind = idx(i,1):idx(i,2);
        w = HTA.ad6.w(ind,bin);
        if rem(length(w),2) == 1
            n = length(w)-1; %force inputs to be even
            w = w(1:end-1,:);
        else
            n = length(w);
        end
        w = detrend(w);
        nfft = n/2;
        window = hanning(n);
        noverlap = nfft/2;
        [Sz1,f] = pwelch(w(:,1),window,noverlap,nfft,fs); %compute PSD
        [Sz2,~] = pwelch(w(:,2),window,noverlap,nfft,fs);
        [Sz3,~] = pwelch(w(:,3),window,noverlap,nfft,fs);
        fr = find(f >= minf & f <= maxf); %set frequency cutoff
        Sz1 = Sz1(fr);
        Sz2 = Sz2(fr);
        Sz3 = Sz3(fr);
        HTA.ad6.Sz1(:,i) = downsample(Sz1,10);
        HTA.ad6.Sz2(:,i) = downsample(Sz2,10);
        HTA.ad6.Sz3(:,i) = downsample(Sz3,10);
        f = f(fr);
        HTA.ad6.f(:,i) = downsample(f,10);
    end
    %1 2 and 3 correspond to bins 1 8 and 16
    HTA.ad6.depths = HTA.ad6.rbdepth(bin);
    
    %AD5117
    sid  = zeros(1,3);
    for i = 1:3
       sdx = find(abs(HTA.ad7.time-samples(i))<=sr); 
       if length(sdx) > 1
           sid(i) = sdx(1);
       else
           sid(i) = sdx;
       end
    end
    idx = [sid(1) sid(1)+avt; sid(2) sid(2)+avt; sid(3) sid(3)+avt];
    for i = 1:length(idx)
        ind = idx(i,1):idx(i,2);
        w = HTA.ad7.w(ind,bin);
        if rem(length(w),2) == 1
            n = length(w)-1; %force inputs to be even
            w = w(1:end-1,:);
        else
            n = length(w);
        end
        w = detrend(w);
        nfft = n/2;
        window = hanning(n);
        noverlap = nfft/2;
        [Sz1,f] = pwelch(w(:,1),window,noverlap,nfft,fs); %compute PSD
        [Sz2,~] = pwelch(w(:,2),window,noverlap,nfft,fs);
        [Sz3,~] = pwelch(w(:,3),window,noverlap,nfft,fs);
        fr = find(f >= minf & f <= maxf); %set frequency cutoff
        Sz1 = Sz1(fr);
        Sz2 = Sz2(fr);
        Sz3 = Sz3(fr);
        HTA.ad7.Sz1(:,i) = downsample(Sz1,10);
        HTA.ad7.Sz2(:,i) = downsample(Sz2,10);
        HTA.ad7.Sz3(:,i) = downsample(Sz3,10);
        f = f(fr);
        HTA.ad7.f(:,i) = downsample(f,10);
    end
    HTA.ad7.depths = HTA.ad7.rbdepth(bin);
    disp('Spectra run. Saving files...')
    save([savedatdir 'VTAvelocities'],'VTA','-v7.3') %note the VTA vecpro file is ~10x larger than the HTA file
    save([savedatdir 'HTAvelocities'],'HTA','-v7.3')
end
clearvars -except HTA VTA c

%Ugh... longest script ever
%It's plot time!!!
f3 = figure(3);
set(f3,'PaperOrientation','portrait',...
    'position',[400 100   800 1000]);
set(gcf,'color','w','PaperPositionMode','auto')
%plot Vertical array of HTA: vpro, vector, vector
s(1) = subplot(311);
[~,n] = size(HTA.topvec.Sz);
for i = 1:n
    p(i) = loglog(HTA.topvec.f(:,i),HTA.topvec.Sz(:,i),...
        'Color',c(i,:),'LineWidth',1.5);hold on
end
int = 1E-2;
xs = linspace(1,9,length(HTA.topvec.Sz));
ys = int.*(xs.^(-3));
p(i+1) = loglog(xs,ys,'Color','k','LineWidth',1.5);
text(3,1E-3,'\itf^-^3')
hold off
g(1) = gca;

s(2) = subplot(312);
[~,n] = size(HTA.midvec.Sz);
for i = 1:n
    p(i) = loglog(HTA.midvec.f(:,i),HTA.midvec.Sz(:,i),...
        'Color',c(i,:),'LineWidth',1.5);hold on
end
int = 1E-2;
xs = linspace(1,9,length(HTA.topvec.Sz));
ys = int.*(xs.^(-3));
p(i+1) = loglog(xs,ys,'Color','k','LineWidth',1.5);
text(3,1E-3,'\itf^-^3')
hold off
g(2) = gca;

s(3) = subplot(313);
[~,n] = size(HTA.vpro.Sz);
for i = 1:n
    p(i) = loglog(HTA.vpro.f(:,i),HTA.vpro.Sz(:,i),...
        'Color',c(i,:),'LineWidth',1.5);hold on
end
int = 1E-2;
xs = linspace(1,9,length(HTA.topvec.Sz));
ys = int.*(xs.^(-5/3));
p(i+1) = loglog(xs,ys,'Color','k','LineWidth',1.5);
text(3,1E-2,'\itf^-^5^/^3')
hold off
g(3) = gca;

set(g(1:3),'Xlim',[1E-1 1E1],'Ylim',[1E-8 1E0],'LineWidth',1.5,'FontSize',12)
% set(g(1:2),'XTickLabel',[])
title(g(1),'\bf\ith = 0.930m','FontSize',12)
title(g(2),'\bf\ith = 0.503m','FontSize',12)
title(g(3),'\bf\ith = 0.185m','FontSize',12)
ylabel(g(2),'\bf\itS(m^-^2s^-^1)','FontSize',12)
xlabel(g(3),'\bf\itFrequency (Hz)','FontSize',12)


f4 = figure(4);
set(f4,'PaperOrientation','portrait',...
    'position',[400 100   800 1000]);
set(gcf,'color','w','PaperPositionMode','auto')
%plot array of HTA aquadopps: bins 1, 8 and 16
s(1) = subplot(321);
[~,n] = size(HTA.ad6.Sz1);
for i = 1:n
    p(i) = loglog(HTA.ad6.f(:,i),HTA.ad6.Sz1(:,i),...
        'Color',c(i,:),'LineWidth',1.5);hold on
end
int = 1E-2;
xs = linspace(1,9,length(HTA.ad6.Sz1));
ys = int.*(xs.^(-5/3));
p(i+1) = loglog(xs,ys,'Color','k','LineWidth',1.5);
text(1,1E-3,'\itf^-^5^/^3')
hold off
g(1) = gca;

s(2) = subplot(323);
[~,n] = size(HTA.ad6.Sz2);
for i = 1:n
    p(i) = loglog(HTA.ad6.f(:,i),HTA.ad6.Sz2(:,i),...
        'Color',c(i,:),'LineWidth',1.5);hold on
end
int = 1E-2;
xs = linspace(1,9,length(HTA.ad6.Sz2));
ys = int.*(xs.^(-5/3));
p(i+1) = loglog(xs,ys,'Color','k','LineWidth',1.5);
text(3,1E-2,'\itf^-^5^/^3')
hold off
g(2) = gca;

s(3) = subplot(325);
[~,n] = size(HTA.ad6.Sz1);
for i = 1:n
    p(i) = loglog(HTA.ad6.f(:,i),HTA.ad6.Sz3(:,i),...
        'Color',c(i,:),'LineWidth',1.5);hold on
end
int = 1E-2;
xs = linspace(1,9,length(HTA.ad6.Sz3));
ys = int.*(xs.^(-5/3));
p(i+1) = loglog(xs,ys,'Color','k','LineWidth',1.5);
text(3,1E-2,'\itf^-^5^/^3')
hold off
g(3) = gca;

s(4) = subplot(322);
[~,n] = size(HTA.ad7.Sz1);
for i = 1:n
    p(i) = loglog(HTA.ad7.f(:,i),HTA.ad7.Sz1(:,i),...
        'Color',c(i,:),'LineWidth',1.5);hold on
end
int = 1E-2;
xs = linspace(1,9,length(HTA.ad7.Sz1));
ys = int.*(xs.^(-5/3));
p(i+1) = loglog(xs,ys,'Color','k','LineWidth',1.5);
text(1,1E-3,'\itf^-^5^/^3')
hold off
g(4) = gca;

s(5) = subplot(324);
[~,n] = size(HTA.ad7.Sz2);
for i = 1:n
    p(i) = loglog(HTA.ad7.f(:,i),HTA.ad7.Sz2(:,i),...
        'Color',c(i,:),'LineWidth',1.5);hold on
end
int = 1E-2;
xs = linspace(1,9,length(HTA.ad7.Sz2));
ys = int.*(xs.^(-5/3));
p(i+1) = loglog(xs,ys,'Color','k','LineWidth',1.5);
text(3,1E-2,'\itf^-^5^/^3')
hold off
g(5) = gca;

s(6) = subplot(326);
[~,n] = size(HTA.ad7.Sz3);
for i = 1:n
    p(i) = loglog(HTA.ad7.f(:,i),HTA.ad7.Sz3(:,i),...
        'Color',c(i,:),'LineWidth',1.5);hold on
end
int = 1E-2;
xs = linspace(1,9,length(HTA.ad7.Sz3));
ys = int.*(xs.^(-5/3));
p(i+1) = loglog(xs,ys,'Color','k','LineWidth',1.5);
text(3,1E-2,'\itf^-^5^/^3')
hold off
g(6) = gca;

set(g(1:6),'Xlim',[1E-1 1E1],'Ylim',[1E-8 1E0],'LineWidth',1.5,'FontSize',12)
% set([g(1:2) g(4:5)],'XTickLabel',[])
% set(g(4:6),'YTickLabel',[])
title(g(1),['\bf\itNo Pneumatophores - h = ' num2str(HTA.ad6.depths(1)) 'm'],'FontSize',12)
title(g(2),['\bf\ith = ' num2str(HTA.ad6.depths(2)) 'm'],'FontSize',12)
title(g(3),['\bf\ith = ' num2str(HTA.ad6.depths(3)) 'm'],'FontSize',12)
title(g(4),['\bf\itPneumatophores - h = ' num2str(HTA.ad7.depths(1)) 'm'],'FontSize',12)
title(g(5),['\bf\ith = ' num2str(HTA.ad7.depths(2)) 'm'],'FontSize',12)
title(g(6),['\bf\ith = ' num2str(HTA.ad7.depths(3)) 'm'],'FontSize',12)
ylabel(g(2),'\bf\itS(m^-^2s^-^1)','FontSize',12)
xlabel(g(3),'\bf\itFrequency (Hz)')
xlabel(g(6),'\bf\itFrequency (Hz)')

