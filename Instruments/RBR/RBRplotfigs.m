%Load and plot RBR files. Will work for either Pressure or Temperature
%Sensors
clear

dirc = 'C:\Users\bkn5\Projects\Mekong_W2015\Data\RBR\SoloT\F2F2';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
cd(dirc)
filelist = dir('RBR*');
files = {filelist.name};
for i = 1:length(files)
    load(files{i});
    dnames = ['d' num2str(i)];
    snames = ['s' num2str(i)];    
    
    yyyy = 2014; mm = 09; dd = 29;
    hr = 01; mi = 00; sec = 00;
    start = yearday(yyyy,mm,dd+hr/24+mi/(24*60)+sec/(24*60*60));
    
    dd = 29;hr = 23; mi = 25; sec = 00;
    finish = yearday(yyyy,mm,dd+hr/24+mi/(24*60)+sec/(24*60*60));
    
    dstr = [num2str(dd) num2str(mm) num2str(yyyy)];
    ind = find(RBR.yearday >= start & RBR.yearday <= finish);
    time = RBR.yearday(ind)-floor(start);
    time2 = time.*24;
    date.(dnames) = time2;
    start2 = (start-floor(start))*24;
    finish2 = (finish-floor(finish))*24;
    if isfield(RBR,'Temperature')
        data.(snames) = RBR.Temperature(ind);
        c = [111,213,255;
            172,176,246;
            233,139,236;
            250,156,156;
            253,191,52;
            248,162,0;
            235,69,0;
            222,0,0;
            209,0,0;
            195,0,0]/255;%custom colormap 'colormap.org'
    end
    if isfield(RBR,'Pressure')
        data.(snames) = RBR.Sea_pressure(ind);
        c = winter(i);
    end
end
m1 = figure('PaperOrientation','portrait',...
    'position',[100 100   1250   550]);
for i = 1:length(files)    
    dnames = ['d' num2str(i)];
    snames = ['s' num2str(i)];
    p(i) = plot(date.(dnames),data.(snames));
    f1 = gca;
    set(p(i),'LineWidth',1.25,'Color',c(i,:))
    set(f1,...
    'Color',[.9 .9 .9],...
    'YGrid','on',...
    'TickDir','in',...
    'XColor',[.2 .2 .2], ...
    'YColor',[.2 .2 .2], ...
    'LineWidth',1,...
    'TickLength',[.02 .02],...
    'YMinorTick','on',...
    'XLim',[start2 finish2])
    hold on

    if isfield(RBR,'Temperature')
        title('RBR Solo-T Temperature Variance')
        ylabel('Degrees C')
        xlabel(['Hour on ' num2str(dd) num2str(mm) num2str(yyyy)])
        leg = legend('10cm','20cm','30cm','40cm','50cm','60cm','70cm','80cm','90cm','100cm');
        set(leg,'Location','SouthWest')
    end
    if isfield(RBR,'Pressure')
        title('RBR Solo-P Pressure Variance')
        ylabel('dBar')
        xlabel(['Hour on ' num2str(dd) num2str(mm) num2str(yyyy)])
        leg = legend('77680','77681','77682','77683');
        set(leg,'Location','NorthEast')
    end 
end

% set(gcf,'color','w','PaperPositionMode','auto');
% filen = ['RBR' dstr];
% filep = 'c:\Users\bkn5\Projects\Mekong_F2014\Mekong_FallDep\Figures\RBR_T\';
% export_fig([filep filen],'-png','-m1','-r900','-opengl')
