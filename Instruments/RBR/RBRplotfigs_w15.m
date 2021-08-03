%Load and plot RBR files. Will work for either Pressure or Temperature
%Sensors
%This code should be used for the 2015 data due to slight differences in
%the way these data were saved.

clear
filen = 'RBR_CT_DPS2';
dirc = 'C:\Users\bkn5\Projects\Mekong_W2015\Data\RBR\SoloT\DPS2';
dirc2 = 'C:\Users\bkn5\Projects\Mekong_W2015\Data\OdysseyCT\CTs\DPS2';
heights1rbr = num2str([170;300;420;540;660]);
heights2rbr = num2str([230;340;460;580;710]);
heights1cts = num2str([110;1020]);
heights2cts = num2str([200 1010]);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
cd(dirc)
filelist = dir('RBR*');
files = {filelist.name};
%since the instruments were deployed in two groupings, 5 in each grouping,
%first select the first 5:
for i = 1:5
    load(files{i});
    sn = num2str(RBR.Metadata.serial);
    if isfield(RBR,'Temperature')
        serial = ['T' sn];
        stack1.(serial) = RBR.Temperature;
        date1.(serial) = RBR.Datetime;
        if rem(i,5)==0
            start = date1.(serial)(1);
            finish = date1.(serial)(end);
        end
        c1 = [0,203,255;
            218,131,255;
            255,216,60;
            255,128,0;
            255,0,0]/255;%custom colormap 'colormap.org'%custom colormap 'colormap.org'
    end
    if isfield(RBR,'Pressure')
        serial = ['P' sn];
        stack1.(serial) = RBR.Sea_pressure;
        date1.(serial) = RBR.Datetime;
        if rem(i,5)==0
            start = date1.(serial)(1);
            finish = date1.(serial)(end);
        end
        c1 = winter(i);
    end
end
for i = 6:10
    load(files{i});
    sn = num2str(RBR.Metadata.serial);
    if isfield(RBR,'Temperature')
        serial = ['T' sn];
        stack2.(serial) = RBR.Temperature;
        date2.(serial) = RBR.Datetime;
        if rem(i,10)==0
            start = date2.(serial)(1);
            finish = date2.(serial)(end);
        end
        c2 = flipud([0,203,255;
            218,131,255;
            255,216,60;
            255,128,0;
            255,0,0])/255;%custom colormap 'colormap.org'
    end
    if isfield(RBR,'Pressure')
        serial = ['P' sn];
        stack2.(serial) = RBR.Sea_pressure;
        date2.(serial) = RBR.Datetime;
        if rem(i,10)==0
            start = date2.(serial)(1);
            finish = date2.(serial)(end);
        end
        c2 = winter(i);
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%load CT sensors
cd(dirc2)
filelist = dir('*.mat');
files = {filelist.name};

for i = 1:length(files)
    load(files{i})
    sn = num2str(CT.Metadata.serial);
    serial = ['CT' sn];
    cttemp.(serial) = CT.Temp;
    ctsalt.(serial) = CT.Salt;
    cttime.(serial) = CT.Datetime;
end
c3 = [190,0,255;
0,0,0]/255;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
m1 = figure('PaperOrientation','portrait',...
    'position',[200 100   1500   900]);
ha = tight_subplot(2,2,0.075,0.05,0.05);
axes(ha(1))
field = fieldnames(date1);
% heights2rbr = flipud(heights2rbr);
mm = repmat(' mm',5,1);
p = zeros(10,1);

for i = 1:5
    p(i) = plot(date1.(field{i}),stack1.(field{i}));
    f1 = gca;
    set(p(i),'LineWidth',1.25,'Color',c1(i,:))
    set(f1,...
        'Color',[.9 .9 .9],...
        'YGrid','on',...
        'TickDir','in',...
        'XColor',[.2 .2 .2], ...
        'YColor',[.2 .2 .2], ...
        'LineWidth',1,...
        'TickLength',[.02 .02],...
        'YMinorTick','on',...
        'XLim',[start finish])
    datetickzoom('x','dd HH:MM:SS','keepticks')
    hold on
    
    title('\bf\itRBR Solo-T Loggers 13-17 - Temperature Variance')
    ylabel('\bf\itDegrees C')
    xlabel('\bf\itdd HH:MM:SS')
    leg = legend([heights1rbr mm]);
    set(leg,'Location','NorthWest')

end

axes(ha(2))
field = fieldnames(date2);
for i = [5 4 3 2 1]
    p(i) = plot(date2.(field{i}),stack2.(field{i}));
    f1 = gca;
    set(p(i),'LineWidth',1.25,'Color',c2(i,:))
    set(f1,...
        'Color',[.9 .9 .9],...
        'YGrid','on',...
        'TickDir','in',...
        'XColor',[.2 .2 .2], ...
        'YColor',[.2 .2 .2], ...
        'LineWidth',1,...
        'TickLength',[.02 .02],...
        'YMinorTick','on',...
        'XLim',[start finish])
    datetickzoom('x','dd HH:MM:SS','keepticks')
    hold on

    title('\bf\itRBR Solo-T Loggers 22-18 - Temperature Variance')
    ylabel('\bf\itDegrees C')
    xlabel('\bf\itdd HH:MM:SS')
    leg = legend([heights2rbr mm]);
    set(leg,'Location','NorthWest')

end

axes(ha(3))
field = fieldnames(cttime);
p = zeros(4,1);
mm = repmat(' mm',2,1);
for i = [3 2]
    if i == 3
        k = 1;
    else
        k = 2;
    end
    p(i) = plot(cttime.(field{i}),ctsalt.(field{i}));set(p(i),'LineWidth',2,'Color',c3(k,:))
    f1 = gca;
    set(f1,...
        'Color',[.9 .9 .9],...
        'YGrid','on',...
        'TickDir','in',...
        'XColor',[.2 .2 .2], ...
        'YColor',[.2 .2 .2], ...
        'LineWidth',1,...
        'TickLength',[.02 .02],...
        'YMinorTick','on',...
        'YLim',[0 35],...
        'XLim',[start finish])
datetickzoom('x','dd HH:MM:SS','keepticks')
hold on

title('\bf\itOdyssey CT Loggers 4055 and 2829 - Salinity')
ylabel('\bf\itmS/cm')
xlabel('\bf\itdd HH:MM:SS')
leg = legend([heights1cts mm]);
set(leg,'Location','NorthWest')
end

axes(ha(4))
field = fieldnames(cttime);
ax = zeros(4,1);p = zeros(4,1);g = zeros(4,1);
mm = repmat(' mm',2,1);
for i = [1 4]
    if i == 4
        k = 2;
    else
        k = 1;
    end
    p(i) = plot(cttime.(field{i}),ctsalt.(field{i}));set(p(i),'LineWidth',2,'Color',c3(k,:))
    f1 = gca;
    set(f1,...
        'Color',[.9 .9 .9],...
        'YGrid','on',...
        'TickDir','in',...
        'XColor',[.2 .2 .2], ...
        'YColor',[.2 .2 .2], ...
        'LineWidth',1,...
        'TickLength',[.02 .02],...
        'YMinorTick','on',...
        'YLim',[0 35],...
        'XLim',[start finish])
datetickzoom('x','dd HH:MM:SS','keepticks')
hold on

title('\bf\itOdyssey CT Loggers 4056 and 2828 - Salinity')
ylabel('\bf\itmS/cm')
xlabel('\bf\itdd HH:MM:SS')
leg = legend([heights1cts mm]);
set(leg,'Location','NorthWest')
end

linkaxes([ha(1) ha(2) ha(3) ha(4)],'x')
if ~exist([filen '.png'],'file')
    set(gcf,'color','w','PaperPositionMode','auto');
    filep = 'c:\Users\bkn5\Projects\Mekong_W2015\Figures\Loggers\RBR_CTs\';
    export_fig([filep filen],'-png','-m1','-r900','-opengl')
end
