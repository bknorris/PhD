%Before working on the bulk of the VecPro data processing, we need to plot
%the time-averaged raw data to look for interesting trends. This script
%loads multiple VecPro files into a single structure, A, runs a temporal average,
%then plots the results. 

clear

%navigate to a folder containing multiple vecpro files
fileList = dir('*VP1*.mat'); %specify which instrument to find based on the
%wildcard operator, ex. *VP1*.mat
fileList = {fileList.name};
%get the data
disp('Reading in the data from multiple files')
for ifile = 1:length(fileList),
    A(ifile) = load(fileList{ifile});
    names{ifile,:} = ['V' fileList{ifile}(28:29)];
end

dfn = fieldnames(A(1).Data); %get fieldnames from first structure element
samprate = A(1).Config.sampleRate;
disp('Averaging velocities in time')
STAT = struct();
STAT.rangebins = A(1).Data.Profiles_Range;
avt = 60*10*samprate; %time to average over in seconds, multiplied by the sample rate
for i = 1:length(A)
    for j = 3:6 %velocity indexes
        dat = A(i).Data.(dfn{j});
        [~,m] = size(dat);
        avdat = zeros(1,m);
        var = regexprep(dfn{j},'Profiles_','');
        for k = 1:m
            avdat(k) = sum(dat(:,k))/avt;
        end
        STAT.(names{i}).Timeav.(var) = avdat;
    end
end

%do some other stuff
% clearvars A dat avdat
%plotting
disp('Plotting time-averaged velocity profiles')
c = jet(length(names));
f1 = figure(1);
set(f1,'PaperOrientation','portrait',...
    'position',[200 200   1250   550]);
ax(1) = subplot(141);
b = zeros(1,length(names));
for i = 1:length(names)
    b(i) = plot(STAT.rangebins,STAT.(names{i}).Timeav.VelBeam1);
    view(-90,90) 
    hold on
    set(b(i),'Color',c(i,:),'LineWidth',2)
end
xlabel('Distance from Transducer (m)')
ylabel('Velocity (m/s)')
title('Beam 1 Velocity Profile')
ax(2) = subplot(142);
b = zeros(1,length(names));
for i = 1:length(names)
    b(i) = plot(STAT.rangebins,STAT.(names{i}).Timeav.VelBeam2);
    view(-90,90) 
    hold on
    set(b(i),'Color',c(i,:),'LineWidth',2)
end
xlabel('Distance from Transducer (m)')
ylabel('Velocity (m/s)')
title('Beam 2 Velocity Profile')
ax(3) = subplot(143);
b = zeros(1,length(names));
for i = 1:length(names)
    b(i) = plot(STAT.rangebins,STAT.(names{i}).Timeav.VelBeam3);
    view(-90,90) 
    hold on
    set(b(i),'Color',c(i,:),'LineWidth',2)
end
xlabel('Distance from Transducer (m)')
ylabel('Velocity (m/s)')
title('Beam 3 Velocity Profile')
ax(4) = subplot(144);
b = zeros(1,length(names));
for i = 1:length(names)
    b(i) = plot(STAT.rangebins,STAT.(names{i}).Timeav.VelBeam4);
    view(-90,90) 
    hold on
    set(b(i),'Color',c(i,:),'LineWidth',2)
end
xlabel('Distance from Transducer (m)')
ylabel('Velocity (m/s)')
title('Beam 4 Velocity Profile')
leg = legend(names,'location','northeastoutside');
set(ax,...
    'XTick',[0.04:0.005:0.07],...
    'XLim', [0.04 0.07],...
    'XGrid','off',...
    'GridLineStyle','-',...
    'YDir','normal',...
    'TickDir','in',...
    'TickLength',[.01 .01],...
    'XMinorTick','on',...
    'YMinorTick','off',...
    'LineWidth',1.25)

