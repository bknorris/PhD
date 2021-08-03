clear

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  It's time for... AQDP data processing!
%  Move this script into the folder that contains the text-converted data files
%  Let's do this in order to avoid confusion: one file at a time
%
%  Developed by Benjamin K. Norris and Dr. Julia C Mullarney, University of
%  Waikato, New Zealand c. 2014
%
%  Processing Scheme:
%  1. Load data, alternatively plot basic figures
%  2. Correct Pressure signal (atmospheric and P-T relationship)
%  3. Trim to data start/stop times
%  4. Adjust velocities by T-S correction
%  (For LR only)
%  5. Rotate to ENU and Magnetic North
%  6. Save Data
%  (For HR only)
%  5. Unwrap problematic data
%  6. Rotate to ENU and Magnetic North
%  7. Save Data
%
%  Dependencies: aqdp_read.m,correctatmp.m,aqdpbeam2enu.m,
%  mag2truenorth.m,interpolate_bad_correlations.m,unwrap_w_prof1.m,
%  interp_nan.m,my_running_median.m,aqdpstrack.m,padbursts.m,
%  cmgbridge.m,cmgidgaps.m
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% <<< Change these values >>>
metdir = 'D:/Projects/Mekong_W2015/Data/Weather/';
dirname = 'C:/users/bkn5/Projects/Mekong_W2015/Data/Aquadopp/Raw/'; %file local directory
fname = 'ADLR4_9March2015'; %instrument name used for saving processed data files
heading = 13; %compass heading of instrument measured in the field
ornt = 'DOWN'; %instrument orientation
HAB = 604; %instrument Height Above Bottom (in mm)
lat  = 9.491840442; %instrument latitude in dms (dd,mm,ss)
lon = 106.2439476; %instrument longitude in dms (dd,mm,ss)
magdec = 0.27; %magnetic declination to be applied
depdate = datenum(2015,3,7,10,00,00); %date deployed (incl. time)
recdate = datenum(2015,3,8,9,00,00); %date recovered (incl. time)

%% Load the data
cd(dirname)
files = dir([fname '*']);
filelist = {files.name};varname = filelist{1};
newvarname = regexprep(varname,'.a1','');

aqdp = aqdp_read(newvarname); %read the data
save([fname '_ver1'],'aqdp','-v7.3') %save intermediate file; ver1 for "version 1"

%% Extract Start/Stop times from the Metadata Record

if strfind(aqdp.metadata.starttime,'a.m.') > 1
    starttime = regexprep(aqdp.metadata.starttime,'a.m.','AM');
    date = regexp(starttime, '(\d+)','match');
    yyyy = str2num(date{3});mm = str2num(date{2});dd = str2num(date{1});
    hh = str2num(date{4});mi = str2num(date{5});sec = str2num(date{6});
end
if strfind(aqdp.metadata.starttime,'p.m.') > 1
    starttime = regexprep(aqdp.metadata.starttime,'p.m.','PM');
    date = regexp(starttime, '(\d+)','match');
    yyyy = str2num(date{3});mm = str2num(date{2});dd = str2num(date{1});
    hh = str2num(date{4});mi = str2num(date{5});sec = str2num(date{6});
    if hh ~= 12;
        hh = hh+12;
    end
end
aqdpstart = datenum(yyyy,mm,dd,hh,mi,sec);

if strfind(aqdp.metadata.endtime,'a.m.') > 1
    endtime = regexprep(aqdp.metadata.endtime,'a.m.','AM');
    date = regexp(endtime, '(\d+)','match');
    yyyy = str2num(date{3});mm = str2num(date{2});dd = str2num(date{1});
    hh = str2num(date{4});mi = str2num(date{5});sec = str2num(date{6});
end
if strfind(aqdp.metadata.endtime,'p.m.') > 1
    endtime = regexprep(aqdp.metadata.endtime,'p.m.','PM');
    date = regexp(endtime, '(\d+)','match');
    yyyy = str2num(date{3});mm = str2num(date{2});dd = str2num(date{1});
    hh = str2num(date{4});mi = str2num(date{5});sec = str2num(date{6});
    if hh ~= 12;
        hh = hh+12;
    end
end
aqdpend = datenum(yyyy,mm,dd,hh,mi,sec);

aqdp.metadata.lat = lat;
aqdp.metadata.lon = lon;
aqdp.metadata.HAB = HAB;

%% Load in MET data, adjust Aquadopp Pressure signal & P-T Adjustment
if 0 %turn flag off (set to 0) if no met data for the record exist
    load([metdir 'Mekong2015_metstation.mat'])
    met.pressure = met.pressure.*0.0001; %convert Pa to dBar
    ind = find(met.datetime >= aqdpstart & met.datetime <= aqdpend);
    met.aqdptime = met.datetime(ind);
    met.aqdppres = met.pressure(ind);
    aqdp.pressure = correctatmp(met.aqdptime,met.aqdppres,aqdp.datenum,aqdp.pressure,aqdp.temperature);
end
%% Clip Start and End Times

%let the user define the start/stop times of the deployment
pause(4)
disp('Clipping datafile to the length of deployment')
int = find(aqdp.datenum >= depdate & aqdp.datenum <= recdate);
snames = fieldnames(aqdp);
if isfield(aqdp,'cor1') %indices of fieldnames differ between LR and HR
    for i = [3:22 25:27]; %clip all fields to start/stop times
        aqdp.(snames{i}) = aqdp.(snames{i})(int,:);
    end
else
    for i = [3:17 20:22]; %clip all fields to start/stop times
        aqdp.(snames{i}) = aqdp.(snames{i})(int,:);
    end
end

%% Apply Sound Speed Correction

t = aqdp.temperature./10;
S = str2num(aqdp.metadata.salinity);
c = 1449.05 + 45.7.*t - 5.21.*t.^2 + 0.23.*t.^3 + (1.333 - 0.126.*t + 0.009.*t.^2).*(S - 35);
%t = T/10 where T = temperature in degrees Celsius, D = 0;
%From:
%A.B. Coppens, Simple equations for the speed of sound in Neptunian waters (1981)
%J. Acoust. Soc. Am. 69(3), pp 862-863
ssnew = mean(c);
ssold = mean(aqdp.sound_spd);
aqdp.vel_b1 = aqdp.vel_b1.*(ssnew/ssold);
aqdp.vel_b2 = aqdp.vel_b2.*(ssnew/ssold);
aqdp.vel_b3 = aqdp.vel_b3.*(ssnew/ssold);
disp('T-S correction for velocities applied')
clear ssnew ssold

%% Sample rate dependencies

if strcmp(aqdp.metadata.samprate,'1 Hz') == 1 %i.e. inst. recorded in low res
    sampmode = 'CONTINUOUS';
else
    sampmode = 'BURST';
end
switch sampmode
    case 'CONTINUOUS'
        % Rotate Velocity Data to ENU then to True North
        switch aqdp.metadata.coordsys
            case 'BEAM'
                %check to see if the internal compass is accurate to the field
                %measurement:
                hcheck = mean(aqdp.heading);
                h1 = hcheck - 3; %+/-3 degrees of inaccuracy
                h2 = hcheck + 3;
                if (heading >= h1 & heading <= h2) == 1
                    disp(['Instrument compass heading ' num2str(floor(hcheck)) ' is accurate to the measured value ' num2str(heading)])
                    heading = aqdp.heading;
                else
                    disp(['Instrument compass heading '  num2str(floor(hcheck)) ' is innacurate, using measured value ' num2str(heading) ' instead'])
                end
                [u,v,w] = aqdpbeam2enu(aqdp.transfm,aqdp.nCells,aqdp.vel_b1,aqdp.vel_b2,aqdp.vel_b3,heading,aqdp.pitch,aqdp.roll,ornt);
            case 'ENU'
                %no rotation is necessary if coordinates are in ENU all
                %ready
                u = aqdp.vel_b1;v = aqdp.vel_b2;w=aqdp.vel_b3;
        end
        %rotate u and v to true north
        [u,v] = mag2truenorth(u,v,magdec);
        %bridge gaps with CMGBRIDGE
        disp('Bridging gaps in velocity time series with CMGBRIDGE')
        disp(['Found ' num2str(cmgidgaps(u)) ' gaps in u time-series'])
        disp(['Found ' num2str(cmgidgaps(v)) ' gaps in v time-series'])
        disp(['Found ' num2str(cmgidgaps(w)) ' gaps in v time-series'])
        nlin = 1E2;
        maxgaps = 1E5;
        u = cmgbridge(u,nlin,maxgaps,maxgaps);
        v = cmgbridge(v,nlin,maxgaps,maxgaps);
        w = cmgbridge(w,nlin,maxgaps,maxgaps);
        aqdp.u = u;aqdp.v = v;aqdp.w = w;
        prompt = 'Run Surface Tracking [y/n]?';
        result = input(prompt,'s');
        if strcmp(result,'y') && strcmp(ornt,'UP') || strcmp(result,'yes') && strcmp(ornt,'UP'); %only run if orientation is UP
            aqdp = aqdpstrack(aqdp,aqdp.metadata.lat);
        elseif strcmp(result,'y') && strcmp(ornt,'DOWN') || strcmp(result,'yes') && strcmp(ornt,'DOWN');
            warning('User opted to run surface tracking on a DOWNWARD facing instrument')
        end
        if strcmp(result,'n') || strcmp(result,'no');
        end
        aqdp.metadata.name = [fname '_f'];
        save(aqdp.metadata.name,'aqdp') %save final file; f for "final"
        disp(['file saved as ' aqdp.metadata.name])
        clearvars -except aqdp heading ornt metdir dirname magdec
        
    case 'BURST'
        % Unrap the data
        disp('Unwrapping the Data: automatic mode')
        % will need to unwrap data burst by burst.
        %parse out a single burst
        spb = aqdp.metadata.spb;
        burstl = (length(aqdp.burst)/spb);
        burstc = ceil(burstl);
        remainder = burstl - burstc; %get the fraction of samples left off at the end
        rind1 = spb*burstc-(spb-1);
        rind2 = spb*burstl;
        %create dummy variables for iterations
        aqdp_burst = struct();
        aqdp.wvel_b1 = [];
        aqdp.wvel_b2 = [];
        aqdp.wvel_b3 = [];
        goodbursts = [];
        badbursts = [];
        
        prompt = 'Skip manual affirmation of auto-unwrapping [y/n]? ';
        result = input(prompt,'s');
        
        if strcmp(result,'y') || strcmp(result,'yes');
           skipunwrapplot = 1;
        end
        if strcmp(result,'n') || strcmp(result,'no');
           skipunwrapplot = 0;
        end
        for i = 1:burstc
            ix1 = spb*i-(spb-1);
            ix2 = spb*i; %this value must always be n+4096 > than ix1
            if i == max(burstc)
                ind = (rind1:rind2);
            else
                ind = (ix1:ix2);
            end
            aqdp_burst.cor1 = aqdp.cor1(ind,:);
            aqdp_burst.cor2 = aqdp.cor2(ind,:);
            aqdp_burst.cor3 = aqdp.cor3(ind,:);
            
            aqdp_burst.vel_b1 = aqdp.vel_b1(ind,:);
            aqdp_burst.vel_b2 = aqdp.vel_b2(ind,:);
            aqdp_burst.vel_b3 = aqdp.vel_b3(ind,:);
            
            aqdp_burst.datenum = aqdp.datenum(ind,:);
            aqdp_burst.pressure = aqdp.pressure(ind,:);
            aqdp_burst.rangebins = aqdp.rangebins;
            
            ccrit=70; %accuracy criteria
            
            % Add in check on percentage of good correlations
            number_goodbeam1=find(aqdp_burst.cor1>=ccrit);
            number_goodbeam2=find(aqdp_burst.cor2>=ccrit);
            number_goodbeam3=find(aqdp_burst.cor3>=ccrit);
            pc_good=sum([length(number_goodbeam1),length(number_goodbeam2),length(number_goodbeam3)])/(aqdp.nCells*aqdp.metadata.spb*3);
            vwrap=(max(aqdp_burst.vel_b1(:))-min(aqdp_burst.vel_b1(:)))*0.5;
            aqdp_burst=interpolate_bad_correlations(aqdp_burst,ccrit);
            
            % Flip only if the instrument is facing upwards
            if strcmp(ornt,'UP') == 1
                aqdp_burst.vel_b1=fliplr(aqdp_burst.vel_b1);
                aqdp_burst.vel_b2=fliplr(aqdp_burst.vel_b2);
                aqdp_burst.vel_b3=fliplr(aqdp_burst.vel_b3);
            end
            
            % Unwrap
            w1=unwrap_w_prof1(aqdp_burst.vel_b1,vwrap);
            w1=fliplr(w1);
            w2=unwrap_w_prof1(aqdp_burst.vel_b2,vwrap);
            w2=fliplr(w2);
            w3=unwrap_w_prof1(aqdp_burst.vel_b3,vwrap);
            w3=fliplr(w3);
            if skipunwrapplot == 0
                % Plot Unwrapped data for inspection
                f1 = figure(9);
                set(f1,'PaperOrientation','portrait',...
                    'position',[100   100   1000   800]);
                a(1) = subplot(311);
                imagesc(aqdp_burst.datenum,aqdp_burst.rangebins, w1')
                hold on
                plot(aqdp_burst.datenum,aqdp_burst.pressure,'k')
                title(['aquadopp burst '  ': VELS UNWRAPPED beam1'])
                c = colorbar;
                caxis([-vwrap vwrap])
                datetick('x','keepticks')
                axis(c);
                ylabel('range (m)')
                a(2) = subplot(312);
                imagesc(aqdp_burst.datenum,aqdp_burst.rangebins,w2')
                hold on
                plot(aqdp_burst.datenum,aqdp_burst.pressure,'k')
                c = colorbar;
                caxis([-vwrap vwrap])
                axis(c);
                title('beam 2')
                datetick('x','keepticks')
                ylabel('range (m)')
                a(3) = subplot(313);
                imagesc(aqdp_burst.datenum,aqdp_burst.rangebins, w3')
                hold on
                plot(aqdp_burst.datenum,aqdp_burst.pressure,'k')
                title('beam 3')
                c = colorbar;
                caxis([-vwrap vwrap])
                axis(c);
                ylabel('range (m)')
                set(a,'YDir','Normal','YLim',[0 0.8])
                datetick('x','keepticks')
                prompt = 'Is the unwrapping satisfactory [y/n]? ';
                result = input(prompt,'s');
                
                if strcmp(result,'y') || strcmp(result,'yes');
                    disp(['Burst ' num2str(i) ' of ' num2str(burstc) ' Saved'])
                    aqdp.wvel_b1(ind,:) = w1;
                    aqdp.wvel_b2(ind,:) = w2;
                    aqdp.wvel_b3(ind,:) = w3;
                    goodbursts(:,i) = i;
                end
                if strcmp(result,'n') || strcmp(result,'no');
                    disp(['Burst ' num2str(i) ' of ' num2str(burstc) ' Not Saved'])
                    badbursts(:,i) = i;
                end
                close(9)
            else
                disp(['Burst ' num2str(i) ' of ' num2str(burstc) ' Saved'])
                aqdp.wvel_b1(ind,:) = w1;
                aqdp.wvel_b2(ind,:) = w2;
                aqdp.wvel_b3(ind,:) = w3;
                goodbursts(:,i) = i;
            end
            clear ix1 ix2 ind int w1 w2 w3 pc_good_sum number_goodbeam1 number_goodbeam2 number_goodbeam3
            
        end
        goodbursts(goodbursts==0) = [];
        badbursts(badbursts==0) = [];
        clear aqdp_burst w1 w2 w3 pc_good_sum number_goodbeam1 number_goodbeam2 number_goodbeam3
        disp(['Bursts Saved: ' num2str(goodbursts)])
        disp(['Bursts Not saved: ' num2str(badbursts)])
        
        if length(badbursts) >= 1
            prompt = 'Do you wish to manually edit the missing burst data [y/n]?';
            result = input(prompt,'s');
            
            if strcmp(result,'y') || strcmp(result,'yes');
                prompt2 = 'Please select a manual edit mode: ';
                disp('Option 1: Manually Edit with ginput [1]')
                disp('Option 2: Selectively NaN out wrapped data [2]')
                result2 = input(prompt2);
                if result2 == 1;
                    % Manually edit the burst data. Will need to break the
                    % bursts up by beam, then by subsections of each beam for
                    % high enough resolution to unwrap the data
                    for i = 1:length(badbursts)
                        yy = badbursts(i);
                        ix1 = spb*yy-(spb-1);
                        ix2 = spb*yy;
                        ind = (ix1:ix2);
                        disp(['Burst ' num2str(yy)])
                        for ii = 1:8
                            jx1 = ii*(spb/8)-(spb/8)+1;
                            jx2 = ii*(spb/8);
                            ind2 = (jx1:jx2);
                            
                            f1 = figure(10);
                            set(f1,'PaperOrientation','portrait',...
                                'position',[100   100   1000   800]);
                            imagesc(aqdp.vel_b1(ind(ind2),:)')
                            title(['aquadopp burst ' num2str(yy) ' Beam 1, Part ' num2str(ii) ' of 8: Manually Unwrap Vels'])
                            c = colorbar;
                            caxis([-vwrap vwrap])
                            datetick('x','keepticks')
                            axis(c);
                            axis xy
                            ylabel('rangebins')
                            aqdp.wvel_b1(ind(ind2),:) = unwrapman(aqdp.vel_b1(ind(ind2),:),vwrap);
                        end
                        close(10)
                        for ii = 1:8
                            jx1 = ii*(spb/8)-(spb/8)+1;
                            jx2 = ii*(spb/8);
                            ind2 = (jx1:jx2);
                            
                            f1 = figure(11);
                            set(f1,'PaperOrientation','portrait',...
                                'position',[100   100   1000   800]);
                            imagesc(aqdp.vel_b1(ind(ind2),:)')
                            title(['aquadopp burst ' num2str(yy) ' Beam 2, Part ' num2str(ii) ' of 8: Manually Unwrap Vels'])
                            c = colorbar;
                            caxis([-vwrap vwrap])
                            datetick('x','keepticks')
                            axis(c);
                            axis xy
                            ylabel('rangebins')
                            aqdp.wvel_b2(ind(ind2),:) = unwrapman(aqdp.vel_b2(ind(ind2),:),vwrap);
                        end
                        close(11)
                        for ii = 1:8
                            jx1 = ii*(spb/8)-(spb/8)+1;
                            jx2 = ii*(spb/8);
                            ind2 = (jx1:jx2);
                            
                            f1 = figure(12);
                            set(f1,'PaperOrientation','portrait',...
                                'position',[100   100   1000   800]);
                            imagesc(aqdp.vel_b1(ind(ind2),:)')
                            title(['aquadopp burst ' num2str(yy) ' Beam 3, Part ' num2str(ii) ' of 8: Manually Unwrap Vels'])
                            c = colorbar;
                            caxis([-vwrap vwrap])
                            datetick('x','keepticks')
                            axis(c);
                            axis xy
                            ylabel('rangebins')
                            aqdp.wvel_b3(ind(ind2),:) = unwrapman(aqdp.vel_b3(ind(ind2),:),vwrap);
                        end
                        close(12)
                    end
                    
                end
                if result2 == 2;
                    %Ask user to pick above which bins to NaN out
                    disp('Program will NaN out bins above the given number')
                    for i = 1:length(badbursts)
                        f1 = figure(13);
                        set(f1,'PaperOrientation','portrait',...
                            'position',[100   100   1000   800]);
                        imagesc(aqdp.vel_b1(ind(ind2),:)')
                        title(['aquadopp burst ' num2str(yy)])
                        c = colorbar;
                        caxis([-vwrap vwrap])
                        axis(c);
                        axis xy
                        ylabel('rangebins')
                        prompt3 = 'Enter a bin number' ;
                        bin = input(prompt3);
                        yy = badbursts(i);
                        ix1 = spb*yy-(spb-1);
                        ix2 = spb*yy;
                        ind = (ix1:ix2);
                        aqdp.wvel_b1(ind,bin:end) = NaN;
                        aqdp.wvel_b2(ind,bin:end) = NaN;
                        aqdp.wvel_b3(ind,bin:end) = NaN;
                        close(13)
                    end
                end
                %If no, fill with NaNs to keep the data size the same as the
                %original unwrapped data
            end
            if strcmp(result,'n') || strcmp(result,'no');
                for i = 1:length(badbursts)
                    yy = badbursts(i);
                    ix1 = spb*yy-(spb-1);
                    ix2 = spb*yy;
                    ind = (ix1:ix2);
                    aqdp.wvel_b1(ind,:) = NaN;
                    aqdp.wvel_b2(ind,:) = NaN;
                    aqdp.wvel_b3(ind,:) = NaN;
                    disp(['Burst ' num2str(yy) ' filled with NaNs'])
                end
            end
        end
        
        disp('Data unwrapping complete')
        save([fname '_ver2'],'aqdp','-v7.3') %save second intermediate file; ver2 for "version 2"
        disp(['file saved as ' fname '_ver2.mat'])
        
        % Rotate Velocity Data to ENU then to True North
        switch aqdp.metadata.coordsys
            case 'BEAM'
                %check to see if the internal compass is accurate to the field
                %measurement:
                hcheck = mean(aqdp.heading);
                h1 = hcheck - 3; %+/-3 degrees of inaccuracy
                h2 = hcheck + 3;
                if (heading >= h1 & heading <= h2) == 1
                    disp(['Instrument compass heading ' num2str(floor(hcheck)) ' is accurate to the measurement ' num2str(heading)])
                    heading = mean(aqdp.heading);
                else
                    disp(['Instrument compass heading '  num2str(floor(hcheck)) ' is innacurate, using measured value ' num2str(heading) ' instead'])  
                end
                [u,v,w] = aqdpbeam2enu(aqdp.transfm,aqdp.nCells,aqdp.wvel_b1,aqdp.wvel_b2,aqdp.wvel_b3,heading,aqdp.pitch,aqdp.roll,ornt);
            case 'ENU'
                %no rotation is necessary if coordinates are in ENU all
                %ready
                u = aqdp.wvel_b1;v = aqdp.wvel_b2;w=aqdp.wvel_b3;
        end
        %rotate u and v to true north
        [u,v] = mag2truenorth(u,v,magdec);
        aqdp.u = u;aqdp.v = v;aqdp.w = w;
        prompt = 'Run Surface Tracking [y/n]?';
        result = input(prompt,'s');
        if strcmp(result,'y') && strcmp(ornt,'UP') || strcmp(result,'yes') && strcmp(ornt,'UP'); %only run if orientation is UP
            aqdp = aqdpstrack(aqdp,aqdp.metadata.lat);
        elseif strcmp(result,'y') && strcmp(ornt,'DOWN') || strcmp(result,'yes') && strcmp(ornt,'DOWN');
            warning('User opted to run surface tracking on a DOWNWARD facing instrument')
        end
        if strcmp(result,'n') || strcmp(result,'no');
        end
        %delete redundant data to reduce file size
        fields = {'vel_b1','vel_b2','vel_b3','wvel_b1','wvel_b2','wvel_b3'};
        aqdp = rmfield(aqdp,fields);
        aqdp.metadata.name = [fname '_f'];
        save(aqdp.metadata.name,'aqdp','-v7.3') %save final file; f for "final"
        disp(['file saved as ' aqdp.metadata.name])
        aqdp = padbursts(aqdp);
        
        %bridge gaps with CMGBRIDGE
        disp('Bridging gaps in velocity time series with CMGBRIDGE')
        disp(['Found ' num2str(cmgidgaps(aqdp.u)) ' gaps in u time-series'])
        disp(['Found ' num2str(cmgidgaps(aqdp.v)) ' gaps in v time-series'])
        disp(['Found ' num2str(cmgidgaps(aqdp.w)) ' gaps in v time-series'])
        nlin = 1E2;
        maxgaps = 1E5;
        aqdp.u = cmgbridge(aqdp.u,nlin,maxgaps,maxgaps);
        aqdp.v = cmgbridge(aqdp.v,nlin,maxgaps,maxgaps);
        aqdp.w = cmgbridge(aqdp.w,nlin,maxgaps,maxgaps);
        disp('Gaps filled')
        
        %Save final file
        fname = [aqdp.metadata.name '_pad'];
        save(fname,'aqdp','-v7.3')
        disp(['file ' fname ' saved'])
        clearvars -except aqdp heading ornt metdir dirname magdec
end