function [Data, Config] = VPQC(Data,Config,brgaps,despike,pflag)
%function to provide preliminary quality control for Vectrino II Profiler
%data. This function analyzes Vectrino data bin by bin first for low
%correlations (<40%) and SNR values (<15dB). Gaps in the data are linearly
%or spectrally filled as NaNs if brgaps is selected.
%Spikes in velocity are removed using the signal processing technique
%developed by Goring & Nikora, 2002, modified by Nobuhito Mori (2005).
%For reference:
%http://www.mathworks.com/matlabcentral/fileexchange/15361-despiking
%
%
% Inputs:
%         Data:    The raw data structure from any .mat-converted Vectrino
%                  file
%         Config:  The raw configuration structure from any .mat-converted
%                  Vectrino file
%         brgaps:  Use a value of 1 for run cmgbridge to linearly/spectrally
%                  bridge gaps in the data record. Use a value of 2 to
%                  alternatively fill gaps with Zeros if this is required by
%                  data processing.
%         despike: Use a value of 1 to run despike routine. Set [] if
%                  undesired.
%         pflag:   [Optional] Boolean flag to plot results, 0 off 1 on
% Outputs:
%         Data:    Modified data structure containing QC'd velocity data
%
%         Config:  A flag is added to the Config file to show QC was applied
%
% NOTE: This script should be run on Vectrino data in BEAM coordinates
% only.
%
% This script was written by Benjamin K Norris, 2015
% University of Waikato, New Zealand
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('Running Vectrino Quality Control')

if nargin < 3
    disp('VPQC requires three inputs')
    help(mfilename);
    return
end

if Config.coordSystem == 1 %1 = XYZ 2 = BEAM
    disp('WARNING: Vectrino velocities are in XYZ coordinates')
    disp('Please convert velocities to BEAM coordinates before running this script')
    return
end

[m,n] = size(Data.Profiles_CorBeam1);
beam1 = zeros(m,n);beam2 = zeros(m,n);beam3 = zeros(m,n);beam4 = zeros(m,n);
V1 = zeros(m,n);V2 = zeros(m,n);V3 = zeros(m,n);V4 = zeros(m,n);

for i = 1:n
    disp(['Bin ' num2str(i)])
    B1 = Data.Profiles_VelBeam1(:,i);B2 = Data.Profiles_VelBeam2(:,i);
    B3 = Data.Profiles_VelBeam3(:,i);B4 = Data.Profiles_VelBeam4(:,i);
    
    %filter low correlations
    ccutoff = 40;
    ind = find(Data.Profiles_CorBeam1(:,i)<ccutoff | Data.Profiles_CorBeam2(:,i)<ccutoff | Data.Profiles_CorBeam3(:,i)<ccutoff | Data.Profiles_CorBeam4(:,i)<ccutoff);
    B1(ind) = NaN;
    B2(ind) = NaN;
    B3(ind) = NaN;
    B4(ind) = NaN;
    
    %filter low SNR values
    acutoff = 15;
    ind = find(Data.Profiles_SNRBeam1(:,i)<acutoff | Data.Profiles_SNRBeam2(:,i)<acutoff | Data.Profiles_SNRBeam3(:,i)<acutoff | Data.Profiles_SNRBeam4(:,i)<acutoff);
    B1(ind) = NaN;
    B2(ind) = NaN;
    B3(ind) = NaN;
    B4(ind) = NaN;
    beam1(:,i) = B1;beam2(:,i) = B1;beam3(:,i) = B3;beam4(:,i) = B4;
end
clear B1 B2 B3 B4
B1 = zeros(m,n);B2 = zeros(m,n);B3 = zeros(m,n);B4 = zeros(m,n);
%report number of bad values found by script
invalid = find(isnan(beam1) | isnan(beam2) | isnan(beam3) | isnan(beam4));
disp(['Eliminated ' num2str(length(invalid)) ' bad values in velocity data']);
disp('Velocity measurements corresponding to correlations <40% have been removed')
disp('Velocity measurements corresponding to SNR values <15dB have been removed')
tic
if brgaps == 1
    disp('Bridging gaps in velocity time series')
    samprate = 50; %in Hz (samples/second)
    nlin = samprate*30; %# of samples in 30sec
    nmaxbr = nlin*8; %# of samples in 4 min
    burst = 10*samprate*60; %# of samples in 10 minutes
    gap1 = cmgidgaps(beam1);gap2 = cmgidgaps(beam2);
    gap3 = cmgidgaps(beam3);gap4 = cmgidgaps(beam4);
    if any([gap1 gap2 gap3 gap4] > 1)
        %crop beams to the same size;
        beam1 = beam1(1:m,1:n);
        beam2 = beam2(1:m,1:n);
        beam3 = beam3(1:m,1:n);
        beam4 = beam4(1:m,1:n);
        id = [1:burst:m m];
        for j = 1:n
            disp(['Bin ' num2str(j)])
            for i = 1:length(id)-1
                B1(id(i):id(i+1),j) = cmgbridge(beam1(id(i):id(i+1),j),nlin,nmaxbr,burst);
                B2(id(i):id(i+1),j) = cmgbridge(beam2(id(i):id(i+1),j),nlin,nmaxbr,burst);
                B3(id(i):id(i+1),j) = cmgbridge(beam3(id(i):id(i+1),j),nlin,nmaxbr,burst);
                B4(id(i):id(i+1),j) = cmgbridge(beam4(id(i):id(i+1),j),nlin,nmaxbr,burst);
            end
        end
    end
    %sometimes cmgbridge creates really large values (e.g. 1E14), find
    %these and set them equal to zero
    B1(B1>10) = 0;B1(B1<-10) = 0;B2(B2>10) = 0;B2(B2<-10) = 0;
    B3(B3>10) = 0;B3(B3<-10) = 0;B4(B4>10) = 0;B4(B4<-10) = 0;
elseif brgaps == 2
    invalid = find(isnan(beam1) | isnan(beam2) | isnan(beam3) | isnan(beam4));
    beam1(invalid) = 0;beam2(invalid) = 0;beam3(invalid) = 0;beam4(invalid) = 0;
    B1 = beam1;B2 = beam2;B3 = beam3;B4 = beam4;
end
clear beam1 beam2 beam3 beam4
disp(['Gaps in data bridged in ' num2str(toc/60) ' minutes'])
%Run Phase Space Thresholding Script
if exist('despike','var')
    disp('Running Phase Space Thresholding Despiking Routine')
    tic
    invalid = find(isnan(B1) | isnan(B2) | isnan(B3) | isnan(B4));
    %PST cannot be run on data with NaNs. Fill with Zeros.
    B1(invalid) = 0;B2(invalid) = 0;B3(invalid) = 0;B4(invalid) = 0;
    for i = 1:n
        disp(['Bin ' num2str(i)])
        [V1(:,i),~] = func_despike_phasespace3d(B1(:,i),[],2);
        [V2(:,i),~] = func_despike_phasespace3d(B2(:,i),[],2);
        [V3(:,i),~] = func_despike_phasespace3d(B3(:,i),[],2);
        [V4(:,i),~] = func_despike_phasespace3d(B4(:,i),[],2);
    end
    disp('Despiking Complete')
    disp(['Velocity data despiked in: ' num2str(toc/60) ' minutes'])
else
    V1 = B1;V2 = B2;V3 = B3;V4 = B4;
end
clear B1 B2 B3 B4
%concatenate velocities into structure for plotting
Dat.V1 = V1;Dat.V2 = V2;Dat.V3 = V3;Dat.V4 = V4;
Org.V1 = Data.Profiles_VelBeam1;Org.V2 = Data.Profiles_VelBeam2; %original
Org.V3 = Data.Profiles_VelBeam3;Org.V4 = Data.Profiles_VelBeam4;
%optional plotting routine
if pflag
    f1 = figure;
    set(f1,'PaperOrientation','portrait',...
        'position',[400 200   1000   800]);
    set(gcf,'color','w','PaperPositionMode','auto')
    sfn = fieldnames(Dat);
    nsamp = 200; %number of bins for histogram
    ui = inputdlg('Enter a bin # to plot (1-35)');
    rb2plot = str2double(ui{1}); %rangebin to plot
    suptitle(['Velocity Distribution, Bin ' num2str(rb2plot) ', n= ' num2str(nsamp)])
    for i = 1:4
        %make xaxis
        xmin = min(Dat.(sfn{i})(:,rb2plot));xmax = max(Dat.(sfn{i})(:,rb2plot));
        xrange = range(Dat.(sfn{i})(:,rb2plot))/nsamp;
        xes1 = xmin:xrange:xmax;
        xmin = min(Org.(sfn{i})(:,rb2plot));xmax = max(Org.(sfn{i})(:,rb2plot));
        xrange = range(Org.(sfn{i})(:,rb2plot))/nsamp;
        xes2 = xmin:xrange:xmax;
        h = hist(Dat.(sfn{i})(:,rb2plot),length(xes2));
        j = hist(Org.(sfn{i})(:,rb2plot),length(xes1));
        %now plot
        ax(i) = subplot(2,2,i);
        p(i) = plot(xes2,j,'-k','linewidth',1);
        hold on
        q(i) = plot(xes1,h,'*r');
        box on
        %         xlabel('\bf\itHz')
        %         ylabel('\bf\itm^2/s')
        title(['\bf\it ' sfn{i}])
        grid on
    end
    leg = legend('\bf\itUnmodified','\bf\itQC Applied','location','northeastoutside');
    xlabel(ax(3),'\bf\itVelocity (m/s)'),xlabel(ax(4),'\bf\itVelocity (m/s)')
    ylabel(ax(1),'\bf\itCounts'),ylabel(ax(3),'\bf\itCounts')
    set(ax(1),'position',[0.1 0.53 0.34 0.35])
    set(ax(2),'position',[0.49 0.53 0.34 0.35])
    set(ax(3),'position',[0.1 0.1 0.34 0.35])
    set(ax(4),'position',[0.49 0.1 0.34 0.35])
    set(leg,'position',[0.85 0.778 0.1 0.1])
    set([ax(1) ax(2) ax(3) ax(4)],...
        'Xlimmode','auto',...
        'GridLineStyle',':',...
        'Xcolor',[0.85 0.85 0.85],...
        'Ycolor',[0.85 0.85 0.85]);
    c_axes = copyobj([ax(1) ax(2) ax(3) ax(4)],gcf);
    set(c_axes, 'color', 'none', 'xcolor', 'k', 'xgrid', 'off', 'ycolor','k', 'ygrid','off');
end
clear Dat Org
%assign variables to Data structure
Data.Profiles_VelBeam1 = V1;
Data.Profiles_VelBeam2 = V2;
Data.Profiles_VelBeam3 = V3;
Data.Profiles_VelBeam4 = V4;
Config.QualityControl = 1; %assign flag to designate that QC has been run
end