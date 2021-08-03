%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Data processing script for the Nortek Vectrino II Profiler. This script
%will run a variety of functions used in post-processing vectrino data. Run 
%this script in the folder that contains the raw data files.
%Select options by turning flags on or off (1 or 0). 

% This script was written by Benjamin K Norris, 2015
% University of Waikato, New Zealand

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear % a good idea
close all
%PREPROCESSING COMMANDS:
concatvp = 1;
crop = 0;
vpqc = 0; %run quality control
cropbottom = 0; %only use if the bottom is visible in the velocity data
%if crop is selected, use these
start = datenum(2014,09,27,13,55,00);
stop = datenum(2014,09,27,15,35,00);
%DATA ANALYSIS COMMANDS:
timeavvp = 0;
velspecvp = 0;
structfunvp = 0;
plotstructfunvp = 0;
plotvp = 0;
savestatfile = 0;
filename = 'VP1_270914';
israwfile = 0; %set this flag to 1 if the file being loaded is a single, raw .mat file
vname = 'VP1'; %vectrino name, either VP1, VP2, or VP3
% savefigdir = 'c:\Users\bkn5\Projects\Mekong_W2015\Figures\Vectrinos\FSS3\VP3\'; 
%fig directory
tic
disp(['Run started at: ' datestr(now)])
%% Load Vectrino Data
if concatvp
    %navigate to a folder containing multiple vecpro files
    fileList = dir(['*' vname '*.mat']); %#ok<*UNRCH> %specify which instrument to find based on the
    %wildcard operator, ex. *VP1*.mat
    fileList = {fileList(3:end).name};
    if exist('filename','var')
        newFile = filename;
    else
        newFile = inputdlg('Enter a file name');
    end
    VPRO = concat_vecpro(fileList,newFile);
else   
    disp(['Loading file ' filename])
    load(filename) %specify concatenated VecPro filename to load
    if israwfile
        VPRO.Data = Data;
        VPRO.Data.Time = Data.Profiles_HostTimeMatlab;
        VPRO.Config = Config;
        clearvars Data Config
    end
end
disp(['Data loaded in: ' num2str(toc/60) ' minutes']),tic
%% define global variables
STAT = struct();
%Start and stop times, user or input data defined:
if exist('start','var') && exist('stop','var')
    crop = 1; %in case you forget to select it
    STAT.DepStart = [datestr(start,'dd-mm-yyyy HH:MM:SS') ' ICT'];
    STAT.DepStop = [datestr(stop,'dd-mm-yyyy HH:MM:SS') ' ICT'];
else
    gmt2ict = (datenum(0,0,0,1,0,0)*7); %ICT = GMT+7; VecPros record in GMT
    STAT.DepStart = [datestr(VPRO.Data.Time(1)+gmt2ict,'dd-mm-yyyy HH:MM:SS') ' ICT'];
    STAT.DepStop = [datestr(VPRO.Data.Time(end)+gmt2ict,'dd-mm-yyyy HH:MM:SS') ' ICT'];
end
STAT.rangebins = VPRO.Data.Profiles_Range;
rb = VPRO.Data.Profiles_Range; %redefine
sr = VPRO.Config.sampleRate;
intv = 10; %averaging window in minutes (must be >= 1 min)
avt = 60*intv*sr;

%% Preliminary processing routines
if crop == 1
    if ~exist('start','var') || ~exist('stop','var')
        warning('Cannot crop timeseries without start and stop times')
        disp('Exiting...')
        return
    end
    disp(['Deployment Start time: ' datestr(start,'dd-mm-yyyy HH:MM:SS')])
    disp(['Deployment Stop time: ' datestr(stop,'dd-mm-yyyy HH:MM:SS')])
    disp('Cropping data to designated start/stop times')
    gmt2ict = (datenum(0,0,0,1,0,0)*7); %ICT = GMT+7; VecPros record in GMT
    if israwfile
        fn = 1:22;
        VPRO.Data.Profiles_HostTimeMatlab = VPRO.Data.Profiles_HostTimeMatlab(1:end-1)+gmt2ict;
        ind = find(VPRO.Data.Profiles_HostTimeMatlab >= start & VPRO.Data.Profiles_HostTimeMatlab <= stop);
    else
        fn = 1:23;
        VPRO.Data.Profiles_HostTimeMatlab = VPRO.Data.Profiles_HostTimeMatlab(1:end-1)+gmt2ict;
        ind = find(VPRO.Data.Profiles_HostTimeMatlab >= start & VPRO.Data.Profiles_HostTimeMatlab <= stop);
    end
    dfn = fieldnames(VPRO.Data);
   for i = fn
       disp(['Cropping field: ' dfn{i}])
       VPRO.Data.(dfn{i}) = VPRO.Data.(dfn{i})(ind,:);
   end
end
if vpqc
    [VPRO.Data,VPRO.Config] = VPQC(VPRO.Data,VPRO.Config,1,2);
end
if cropbottom
    disp('User requests sub bed-level velocity measurements to be removed')
    %first need to refine the bottom dist parameter
    bdist = VPRO.Data.BottomCheck_BottomDistance;
    bdmax = runningmax(bdist,100);
    bdmax = my_running_median(bdmax,500); %despike
    bdmax = smooth(bdmax,100,'sgolay');
    %need to make bdmax the same length as the t-s
    ts = 1:length(VPRO.Data.Profiles_VelBeam1);
    bd = spline(VPRO.Data.BottomCheck_HostTimeMatlab,bdmax,VPRO.Data.Time);
    if length(VPRO.Data.Time) ~= length(VPRO.Data.Profiles_VelBeam1)
        bd = bd(1:length(VPRO.Data.Profiles_VelBeam1));
    end
    dhts = repmat(VPRO.Data.Profiles_Range,length(bd),1);
    for ii = 1:length(bd)
        indx = dhts(ii,:) >= bd(ii);
        VPRO.Data.Profiles_VelBeam1(ii,indx) = NaN;
        VPRO.Data.Profiles_VelBeam2(ii,indx) = NaN;
        VPRO.Data.Profiles_VelBeam3(ii,indx) = NaN;
        VPRO.Data.Profiles_VelBeam4(ii,indx) = NaN;
        count = find(indx == 1);
    end
    %We need to make a mask for the TKE estimates by finding the minimum depth
    %of the bed over the given averaging interval. The no. of lags defines how 
    %many points are subsampled from the dataset for the difference 
    %calculation. We use this number here as it is the maximum depth that the
    %eventual TKE estimates will use.
    lags = 5;
    idx = avt:avt:length(bd);
    idxx = [1 idx];
    k = 1;
    while k < length(idxx)
        idx = idxx(k):idxx(k+1);
        [~,invc] = find(isnan(VPRO.Data.Profiles_VelBeam1(idx,:)));
        if isempty(invc)
            bdmin(k) = 9999; %bunk number in case there are no measurements conflicting with the bottom distance
            warning(['No bottom detected in profile at timestep ' num2str(k) '; (period) ' num2str(k*intv) ' min'])
        else
            bdmin(k) = min(invc);
        end
        k = k+1;
    end
    disp('Sub bed-level velocity measurements have been removed')
end
disp(['Data preprocessed in: ' num2str(toc/60) ' minutes']),tic
%% Generate Time averages for profile plots
if timeavvp
    %First, transform beam -> XYZ
    [VPRO.Data,VPRO.Config] = VPro_coordinateTransform(VPRO.Data,VPRO.Config,'bx');
    tf_flag = 1; %create a flag if the transform was run for future instances
    dfn = fieldnames(VPRO.Data);
    disp('Averaging velocities in time')
    for i = 65:68; %select only velocity fields
        [n,m] = size(VPRO.Data.(dfn{i}));
        idx = avt:avt:n;
        idxx = [1 idx];
        dat = VPRO.Data.(dfn{i});
        %do the averaging
        avdat = zeros(length(idxx)-1,m);
        k = 1;
        while k < length(idxx)
            idx = idxx(k):idxx(k+1);
            for j = 1:m
                avdat(k,j) = (sum(dat(idx,j))/numel(dat(idx,j)));
            end
            k = k+1;
        end
        %sometimes HUGE values are produced, remove these
        if max(max(avdat)) > 5|| min(min(avdat)) < -5 %set limit to greater or less than 5m/s
            avind = find(avdat > 5 | avdat < -5);
            avdat(avind) = 0;
        end
        %name the fields
        var = regexprep(dfn{i},'Profiles_Vel','');
        vari = ['AvV' var];
        STAT.Avs.cmt = 'Each row of data is the average velocity for the given interval';
        STAT.Avs.interval = [num2str(intv) ' minutes'];
        STAT.Avs.(vari) = avdat;
        clearvars dat avdat idx idxx n m var vari
    end
end

%% Calculate Spectra
if velspecvp
    disp('Calculating Velocity Spectra')
    load ab_nu5percent %stored in C:\Users\Documents\MATLAB\
    %First, transform beam -> XYZ
    if ~exist('tf_flag','var')
        [VPRO.Data,VPRO.Config] = VPro_coordinateTransform(VPRO.Data,VPRO.Config,'bx');
        tf_flag = 1; %create a flag if the transform was run for future instances
    end
    if cropbottom %nans out bad values, spectra doesn't play nicely with NaNs
        invalid = find(isnan(VPRO.Data.Profiles_VelX) | isnan(VPRO.Data.Profiles_VelY) |...
            isnan(VPRO.Data.Profiles_VelZ1) | isnan(VPRO.Data.Profiles_VelZ2));
        VPRO.Data.Profiles_VelX(invalid) = 0;
        VPRO.Data.Profiles_VelY(invalid) = 0;
        VPRO.Data.Profiles_VelZ1(invalid) = 0;
        VPRO.Data.Profiles_VelZ2(invalid) = 0;
    end
    dfn = fieldnames(VPRO.Data);
    %define spectra parameters
    nf = sr/2; %nyquist criterium
    nfft = 0.15*avt;
    overlap = 0.7;
    overlpts = floor(overlap*nfft);
    %load in the data separately from processing to concatenate the Z beams
    for i = 65:68
        disp(['Reading field: ' dfn{i}])
        [n,~] = size(VPRO.Data.(dfn{i}));
        idx = avt:avt:n;
        idxx = [1 idx]; %do not include remainders here for this calculation
        dat = VPRO.Data.(dfn{i});
        dat = detrend(dat); %spectra need detrended data
        %preallocate variables
        avdat = zeros(length(idx),1);
        %calculate specta of field
        k = 1;
        while k < length(idxx)
            idx = idxx(k):idxx(k+1);
            %run spectra on depth-averaged velocities
            for j = 1:length(idx)
                avdat(j,:) = mean(dat(idx(j),:));
            end
            [spec(:,k),f(:,k)] = cpsd(avdat,avdat,hanning(nfft,'periodic'),overlpts,nfft,sr); %#ok<*SAGROW>
            k = k+1;
        end
        %calculate degrees of freedom and the confidence intervals
        dof = floor(8*length(avdat)/(3*nfft));
        fprintf('Degrees of freedom: %d\n',dof)
        ciLimit = [0.025,0.975]; % For 95% CIs.
        ab = ab_nu5(dof,:);
        ci_up=dof.*(spec)./ab(1); %Upper
        ci_lw=dof.*(spec)./ab(2); %Lower
        %name the fields
        var = regexprep(dfn{i},'Profiles_','');
        STAT.Spec.cmt = 'Each column of data is the spectra for the given interval';
        STAT.Spec.interval = [num2str(intv) ' minutes'];
        STAT.Spec.(var).Spec = spec;
        STAT.Spec.(var).CI_up = ci_up;
        STAT.Spec.(var).CI_low = ci_lw;
        clearvars dat avdat idx idxx n m var spec
    end
    STAT.Spec.f = f;
    if cropbottom %put the NaNs back for other routines
        VPRO.Data.Profiles_VelX(invalid) = NaN;
        VPRO.Data.Profiles_VelY(invalid) = NaN;
        VPRO.Data.Profiles_VelZ1(invalid) = NaN;
        VPRO.Data.Profiles_VelZ2(invalid) = NaN;
    end
end

%% Run Structure Function - Calculate TKE Dissipation rate
if structfunvp
    %This code borrows the method from Wiles, et al. 2006 "A novel 
    %technique for measuring the rate of turbulent dissipation in the 
    %marine environment"
    disp('Calculating TKE Dissipation rates')
    %structure function cannot be run on beams in XYZ. Transform from XYZ
    %-> BEAM
    if exist('tf_flag','var')
        [VPRO.Data,VPRO.Config] = VPro_coordinateTransform(VPRO.Data,VPRO.Config,'xb');
        fn = 65:68;
    elseif israwfile
        fn = 3:6;
    else
        fn = 4:7;
    end
    if cropbottom %revert bottom measurements back to zero vels for the TKE calculation
        VPRO.Data.Profiles_VelBeam1(invalid) = 0;
        VPRO.Data.Profiles_VelBeam2(invalid) = 0;
        VPRO.Data.Profiles_VelBeam3(invalid) = 0;
        VPRO.Data.Profiles_VelBeam4(invalid) = 0;
    end    
    dfn = fieldnames(VPRO.Data);
    r = (VPRO.Config.cellSizeSelected/1000)/cosd(30); %in m, converts
    %vertical beam distance to along beam distance. This is necessary for
    %the structure function analysis.
    for i = fn
        disp(['Reading field: ' dfn{i}])
        [n,m] = size(VPRO.Data.(dfn{i}));  
        idx = avt:avt:n;
        idxx = [1 idx]; %do not include remainders here for this calculation
        dat = detrend(VPRO.Data.(dfn{i}),'constant');
        %According to Wiles et al., Structure function analysis requires
        %the temporal mean to be removed from the timeseries. The no. of 
        %lags defines how many points are subsampled from the dataset for 
        %the difference calculation.
        if ~exist('lags','var')
            lags = 5;
        end
        %initialize local variables:
        maxbin = m-lags; %TKE estimates are maximally length(maxbin) due to lag #
        E = zeros(maxbin,length(idxx)-1);
        D = zeros(maxbin,lags);
        pf = zeros(maxbin,2);
        pv = zeros(maxbin,lags);
        R = zeros(maxbin,lags);
        rr = zeros(1,lags);
        N = zeros(1,maxbin);
        epsilon = zeros(1,maxbin);
        k = 1;

        %loop in time
        while k < length(idxx) %iteration counter, i.e. num of subdivisions of the entire dataset
            idx = idxx(k):idxx(k+1);
            disp(['Calculating TKE dissipation rates, cycle ' num2str(k) ' of ' num2str(length(idxx)-1)])
            itt = 1; %iteration number
            %loop in depth
            while itt <= maxbin
                d = zeros(length(idx),lags);
                idx2 = itt:itt+lags;
                %iterate velocity differences (vel(z,1) - vel(z,2:lags))^2
                for j = 1:lags
                    d(:,j) = (dat(idx,itt)-(dat(idx,idx2(j+1)))).^2;
                    D(itt,j) = nanmean(d(:,j));
                    rr(:,j) = r*j;
                end
                R(itt,:) = rr.^(2/3); %along beam distance ^ 2/3
                pf(itt,:) = polyfit(R(itt,:),D(itt,:),1); %linear regression of D along r^2/3
                pv(itt,:) = polyval(pf(itt,:),R(itt,:));
                %calculate estimate of epsilon (slope = A =
                %(Cv^2*epsilon^2/3) where Cv^2 = 2.1 (Monin, A.S. and A.M.
                %Yaglom. 1975. Statistical Fluid Mechanics: Mechanics of Turbulence,
                %Volume 2, Dover, New York.)
                A = pf(itt,1);
                N(itt) = pf(itt,2); 
                epsilon(itt) = (A/2.1)^(3/2); %units of W/m^3
                itt = itt+1;
            end
            E(1:length(epsilon),k) = epsilon'; %TKE dissipation rate
            E(abs(imag(E))>0) = NaN; %filter non-real numbers from E estimates, Noise
            E(abs(imag(N'))>0) = NaN;
            if cropbottom
                %mask out below-bottom estimates if cropbottom is on
                for ii = 1:length(bdmin)
                    E(bdmin(ii):end,ii) = NaN;
                end
            end
            %make time string for plots (in minutes)
            Etime(:,k) = intv*k; 
            k = k+1;
        end
        %compute estimates of log10(E) for plotting
        Emean = log10(E);
        %name the fields
        var = regexprep(dfn{i},'Profiles_Vel','');
        STAT.TKE.cmt = 'Each column of data is the TKE dissipation rates for the given interval';
        STAT.TKE.interval = [num2str(intv) ' minutes'];
        STAT.TKE.(var) = E;
    
        %plot results, pick a bin in the middle of the beam to plot D(z,r)
        if plotstructfunvp
            f1 = figure;
            set(f1,'PaperOrientation','portrait',...
                'position',[400 200   1000   600]);
            ax = tight_subplot(1,2,0.1,[0.1 0.05],[0.1 0.05]);
            legfrmt = 'z = %d m';
            counter = 1:3:maxbin;
            c = jet(length(counter));
            axes(ax(1)) %#ok<*LAXES>
            for a = 1:length(counter);
                p1(a) = plot(R(1,:),D(counter(a),:),'Marker','o',...
                    'LineStyle','none','Color',c(a,:));
                hold on
                p2(a) = plot(R(1,:),pv(counter(a),:),'LineWidth',1.5,...
                    'Color',c(a,:));
                hold on
                
                names{a} = sprintf('z = %.3f m',rb(counter(a)));
            end
            %crop output to size(a)
            p2(p2 == 0) = [];
            emptyCells = cellfun('isempty',names); 
            names(emptyCells) = [];
            leg1 = legend(p2,names,'location','northeastoutside');
            ylabel('\bf\itVelocity Difference Squared (m^2s^-^2)','FontSize',14)
            xlabel('\bf\itAlong Beam Distance (m^2^/^3)','FontSize',14)
%             ordmag = 10^floor(log10(max(max(D)))); %calculate order of magnitude to set axis limits
            ulx = max(max(R))+1E-2;uly = max(max(D));
            lly = min(min(D));
            set(gca,...
                'XTick', 1E-2:5E-3:ulx,...
                'TickDir','in',...
                'TickLength',[.01 .01],...
                'XMinorTick','off',...
                'YMinorTick','off',...
                'LineWidth',1.25)
            box on
            axis([1E-2 ulx lly uly])
            axes(ax(2))
            [~,m] = size(Emean);
            c = linspecer(m,'sequential');
            clear names
            for a = 1:m
                p3(a) = plot(Emean(1:maxbin,a),rb(1:maxbin));
                set(p3(a),'LineWidth',1,'Color',c(a,:),'Marker','o')
                hold on
                names{a} = sprintf('%d min to %d min',(intv*a-intv),intv*a);
            end
            llx = min(min(Emean));ulx = max(max(Emean));
            leg2 = legend(p3,names,'location','northeastoutside');
            ylabel('\bf\itDepth (m)')
            xlabel('\bf\itlog_1_0( \epsilon ) (W/kg)','FontSize',14)
            set(gca,'YDir','reverse','LineWidth',1.25)
            set(ax,'FontSize',14,'FontName','Helvetica')
            axis([llx ulx 0.04 0.07])
%             set(ax(1),'position',[0.1 0.1 0.35 0.85])
%             set(ax(2),'position',[0.62 0.1 0.22 0.85])
%             set(leg1,'position',[0.44 0.53 0.1 0.4])
%             set(leg2,'position',[0.87 0.63 0.1 0.3])
            set(gcf,'color','w','PaperPositionMode','auto')
            prompt = 'Save Figure? [y/n] ';
            result = input(prompt,'s');
            if strcmp(result,'y') || strcmp(result,'yes');
                var = regexprep(dfn{i},'Profiles_Vel','');
                fpath = savefigdir;fname = [vname '_sf_p_' var];
                export_fig([fpath fname],'-png','-m1','-r900','-opengl')
                disp(['Figure ' fname '.png saved'])
            end
            if strcmp(result,'n') || strcmp(result,'no');
            end
        end
        %delete extra bulk data
        clearvars D d E dat idx idx2 idxx j k maxbin rr p1 p2 p3 ax Emean pf pv epsilon R A
        clearvars legfrmt itt leg names f1 c b a minl maxl leg1 leg2 ulx uly lly llx ulx Eimag Eind
    end
    %define time vector
    STAT.TKE.Time = Etime;
    if cropbottom %replace NaNs for plot routines as is required. 
        VPRO.Data.Profiles_VelBeam1(invalid) = NaN;
        VPRO.Data.Profiles_VelBeam2(invalid) = NaN;
        VPRO.Data.Profiles_VelBeam3(invalid) = NaN;
        VPRO.Data.Profiles_VelBeam4(invalid) = NaN;
    end   
end


%% Plot
if plotvp
    plotVecProfiles(VPRO.Data,0.25,0.15)
    prompt = 'Save Figure? [y/n] ';
    result = input(prompt,'s');
    if strcmp(result,'y') || strcmp(result,'yes');
        fpath = savefigdir;fname = [vname '_timeseries'];
        export_fig([fpath fname],'-png','-m1','-r900','-opengl')
        disp(['Figure ' fname '.png saved'])
    end
    if strcmp(result,'n') || strcmp(result,'no');
    end
    if timeavvp
        f1 = figure;
        set(f1,'PaperOrientation','portrait',...
            'position',[400 200   800   600]);
        set(gcf,'color','w','PaperPositionMode','auto')
        [n,~] = size(STAT.Avs.AvVX);
        sfn = fieldnames(STAT.Avs);
        c = linspecer(n,'sequential');
        titles = {'X';'Y';'Z1';'Z2'};
        for i = 1:4
            ib = 3:length(sfn);
            ax(i) = subplot(1,4,i);
            for ii = 1:n
                p1(ii) = plot(STAT.Avs.(sfn{ib(i)})(ii,:),rb,'Color',c(ii,:),...
                    'LineWidth',1.5);
                hold on
                xlabel('\bf\itm/s')
                title(['\bf\it ' titles{i}])
                set(gca,'YDir','reverse')
                names{ii} = sprintf('%d min to %d min',(intv*ii-intv),intv*ii);
            end
            xl = xlim;
            axis([xl(1) xl(2) rb(1) rb(end)])
        end
        leg = legend(p1,names,'location','northeastoutside');
        ylabel(ax(1),'\bf\itDepth (m)')
        set([ax(2) ax(3) ax(4)],'YTickLabel',[])
        suptitle(['\bf\itTime-Averaged Velocities, ' vname])
        set([ax(1) ax(2) ax(3) ax(4)],'Xlimmode','auto')
        set(ax(1),'position',[0.1 0.1 0.14 0.777])
        set(ax(2),'position',[0.28 0.1 0.14 0.777])
        set(ax(3),'position',[0.46 0.1 0.14 0.777])
        set(ax(4),'position',[0.64 0.1 0.14 0.777])
        set(leg,'position',[0.84 0.45 0.1 0.1])
        set(gcf,'color','w','PaperPositionMode','auto')
        prompt = 'Save Figure? [y/n] ';
        result = input(prompt,'s');
        if strcmp(result,'y') || strcmp(result,'yes');
            fpath = savefigdir;fname = [vname '_tav_vel'];
            export_fig([fpath fname],'-png','-m1','-r900','-opengl')
            disp(['Figure ' fname '.png saved'])
        end
        if strcmp(result,'n') || strcmp(result,'no');
        end
    end
    if velspecvp
        %plot example spectra
        titles = {'X';'Y';'Z1';'Z2'};
        f1 = figure;
        set(f1,'PaperOrientation','portrait',...
            'position',[400 200   1000   800]);
        set(gcf,'color','w','PaperPositionMode','auto')
        [~,m] = size(STAT.Spec.VelX.Spec);
        c = jet(m);
        sfn = fieldnames(STAT.Spec);
        for i = 1:4
            ax(i) = subplot(2,2,i);
            ib = 3:length(sfn)-1;
            for ii = 1:m
                area(STAT.Spec.f(:,ii),STAT.Spec.(sfn{ib(i)}).CI_up(:,ii),'FaceColor',[0.9 0.9 0.9],'LineStyle','none'),hold on
                area(STAT.Spec.f(:,ii),STAT.Spec.(sfn{ib(i)}).CI_low(:,ii),'FaceColor',[1 1 1],'LineStyle','none'),hold on
                p(ii) = plot(STAT.Spec.f(:,ii),STAT.Spec.(sfn{ib(i)}).Spec(:,ii),'-x','linewidth',1.5);
                set(p(ii),'Color',c(ii,:))
                hold on
                box on
                xlabel('\bf\itHz')
                ylabel('\bf\itm^2/s')
                title(['\bf\it ' titles{i}])
                names{ii} = sprintf('%d min to %d min',(intv*ii-intv),intv*ii);
            end
        end
        leg = legend(p,names,'location','northeastoutside');
        set([ax(1) ax(2)],'Xlim',[0 1],'Ylim',[0 0.01])
        set([ax(3) ax(4)],'Xlim',[0 1],'Ylim',[0 0.01])
        suptitle(['\bf\itXYZ Velocities Spectra, ' vname])
        set(ax(1),'position',[0.1 0.53 0.30 0.35])
        set(ax(2),'position',[0.48 0.53 0.30 0.35])
        set(ax(3),'position',[0.1 0.1 0.30 0.35])
        set(ax(4),'position',[0.48 0.1 0.30 0.35])
        set(leg,'position',[0.84 0.45 0.1 0.1])
        prompt = 'Save Figure? [y/n] ';
        result = input(prompt,'s');
        if strcmp(result,'y') || strcmp(result,'yes');
            fpath = savefigdir;fname = [vname '_spectra'];
            export_fig([fpath fname],'-png','-m1','-r900','-opengl')
            disp(['Figure ' fname '.png saved'])
        end
        if strcmp(result,'n') || strcmp(result,'no');
        end
        %random plots: Spectral Averages
        Xes = STAT.Spec.VelX.Spec+STAT.Spec.VelY.Spec;
        Zes = (STAT.Spec.VelZ1.Spec+STAT.Spec.VelZ2.Spec)./2;
%         intv = 10;
        f1 = figure;
        set(f1,'PaperOrientation','portrait',...
            'position',[400 200   1000   800]);
        set(gcf,'color','w','PaperPositionMode','auto')
        % suptitle('\bf\itS(X+Y) and S(Z1+Z2)/2')
        [~,m] = size(Xes);
        c = jet(m);
        ax(1) = subplot(2,1,1);
        for ii = 1:m
            area(STAT.Spec.f(:,ii),STAT.Spec.VelX.CI_up(:,ii),'FaceColor',[0.9 0.9 0.9],'LineStyle','none'),hold on
            area(STAT.Spec.f(:,ii),STAT.Spec.VelX.CI_low(:,ii),'FaceColor',[1 1 1],'LineStyle','none'),hold on
            p(ii) = plot(STAT.Spec.f(:,ii),Xes(:,ii),'-x','linewidth',1.5);
            set(p(ii),'Color',c(ii,:))
            hold on
            box on
            %     xlabel('\bf\itHz')
            ylabel('\bf\itm^2/s')
            %     title('\bf\itHorizontal Spectra S(X+Y)')
            names{ii} = sprintf('%d min to %d min',(intv*ii-intv),intv*ii);
            title('\bf\itS(X+Y)')
        end
        grid on
        %leg = legend(p,names,'location','northeast');
        ax(2) = subplot(2,1,2);
        for ii = 1:m
            area(STAT.Spec.f(:,ii),STAT.Spec.VelZ1.CI_up(:,ii),'FaceColor',[0.9 0.9 0.9],'LineStyle','none'),hold on
            area(STAT.Spec.f(:,ii),STAT.Spec.VelZ1.CI_low(:,ii),'FaceColor',[1 1 1],'LineStyle','none'),hold on
            p(ii) = plot(STAT.Spec.f(:,ii),Zes(:,ii),'-x','linewidth',1.5);
            set(p(ii),'Color',c(ii,:))
            hold on
            box on
            xlabel('\bf\itHz')
            ylabel('\bf\itm^2/s')
            title('\bf\itS(Z1+Z2)/2')
            %     title('\bf\itAverage Vertical Spectra S(Z1+Z2)/2')
        end
        leg = legend(p,names,'location','northeast');
        grid on
        set([ax(1) ax(2)],...
            'Xlimmode','auto',...
            'GridLineStyle',':')
        set(ax(1),...
            'Xlim',[0 1],...
            'Ylim',[0 0.05],...
            'XTickLabel',[],...
            'YTick',0:0.025:0.05,...
            'position',[0.1 0.52 0.8 0.4])
        set(ax(2),...
            'Xlim',[0 1],...
            'Ylim',[0 0.01],...
            'YTick',0:0.005:0.01,...
            'position',[0.1 0.08 0.8 0.4])
        set(leg,'position',[0.8 0.45 0.1 0.1])
        prompt = 'Save Figure? [y/n] ';
        result = input(prompt,'s');
        if strcmp(result,'y') || strcmp(result,'yes');
            fpath = savefigdir;fname = [vname '_spectraAvgs'];
            export_fig([fpath fname],'-png','-m1','-r900','-opengl')
            disp(['Figure ' fname '.png saved'])
        end
        if strcmp(result,'n') || strcmp(result,'no');
        end
    end
    if structfunvp
        %plot E(t,r)
        %pcolor requests time/depth arrays the same size as the data to be
        %plotted...
        sfn = fieldnames(STAT.TKE);
        for i = 1:4
            ib = 3:length(sfn)-1;
            f1 = figure;
            set(f1,'PaperOrientation','portrait',...
                'position',[400 200   800   800]);
            set(gcf,'color','w','PaperPositionMode','auto')
            [m,n] = size(STAT.TKE.(sfn{ib(i)}));
            data = log10(STAT.TKE.(sfn{ib(i)}));
            mini = round(min(min(data)))-1;
            maxi = round(max(max(data)))+1;
            time = STAT.TKE.Time;
            time = repmat(time,m,1);
            rangebins = repmat(rb(1:m),n,1)';
            p(i) = pcolor(time,rangebins,data);
            shading interp
            colormap jet
            caxis([mini maxi])
            c = colorbar;
            cbh = get(c,'ylabel');
            titleString = '\bf\itlog_1_0(\epsilon) (W/kg)';
            set(cbh ,'String',titleString,'FontSize',14);
            set(c,'YTick',mini:0.5:maxi);
            set(p(i),'linestyle','none')
            set(gca,'layer','top')
            set(gca,...
                'YDir','Reverse',...
                'YTick',0.04:0.005:0.07,...
                'YLim',[0.04 0.07],...
                'YGrid','on',...
                'GridLineStyle','--',...
                'Xlim',[time(1) time(end)],...
                'box','on');
            xlabel('\bf\itTime Elapsed (min)')
            ylabel('\bf\itDepth (m)')
            title(['\bf\itEstimates of Log_1_0(\epsilon) , ' sfn{ib(i)}])
            prompt = 'Save Figure? [y/n] ';
            result = input(prompt,'s');
            if strcmp(result,'y') || strcmp(result,'yes');
                fpath = savefigdir;fname = [vname '_TKE_' sfn{ib(i)}];
                export_fig([fpath fname],'-png','-m1','-r900','-opengl')
                disp(['Figure ' fname '.png saved'])
            end
            if strcmp(result,'n') || strcmp(result,'no');
            end
        end
    end
end
disp(['Data analysis completed in: ' num2str(toc/60) ' minutes'])
if savestatfile
    %save file
    sfname = ['Stat_' filename];
    save(sfname,'STAT','-v7.3')
    disp(['Statistics file saved as ' sfname])
end
clearvars -except STAT VPRO start stop filename concatvp israwfile plotstructfunvp plotvp rb savefigdir savestatfile structfunvp tf_flag timeavvp velspecvp vname