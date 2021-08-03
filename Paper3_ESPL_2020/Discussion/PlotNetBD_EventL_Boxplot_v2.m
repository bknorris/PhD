%Plot the combined event data
clear,close all
load('e:\Mekong_W2015\DataAnalysis\Paper3\Wavelet\CombEventdata_flood_v5.mat');
data.flood = dat;
load('e:\Mekong_W2015\DataAnalysis\Paper3\Wavelet\CombEventdata_ebb_v5.mat');
data.ebb = dat;clear dat
fn = fieldnames(data);
sfdir = 'e:\MemoryStick\GradSchool\DataAnalysis\Paper3\Figures\';
flds = {'mud';'fringe';'forest'};
band = {'wave';'ig'};
cl = [178,24,43;5,55,122]./255;
symb = {'o';'^';'s'};
%% Net Bed Level - Event Length (Acc & Ero)
sp = zeros(3,2);
eb = zeros(2,1);
f1 = figure(1);
set(f1,'PaperOrientation','portrait',...
    'position',[400 100   650   450]);
w = [1 4; 2 5; 3 6];
value = 0; %find values > 1 mm
for j = 1:length(flds)
    
    sp(w(j,1)) = subplot(2,3,w(j,1));
%     wvf = data.flood.(flds{j}).wave.deltbd;
%     wve = data.ebb.(flds{j}).wave.deltbd;
    igf = data.flood.(flds{j}).ig.deltbd;
    ige = data.ebb.(flds{j}).ig.deltbd;
    %Accretion
%     id1 = find(wvf>value);
%     id2 = find(wve>value);
    id3 = find(igf>value);
    id4 = find(ige>value);
    wbd = [igf(id3); ige(id4)];
    grp = [repmat(1,length(igf(id3)),1); repmat(2,length(ige(id4)),1)];
    plot(ones(10,1)*2.5,linspace(-0.05,0.15,10),'-k'),hold on
    bp = boxplot(wbd,grp,'symbol','.',...
        'colors',cl,...
        'boxstyle','outline',...
        'plotstyle','traditional',...
        'jitter',1,...
        'outliersize',11);
    set(bp,'linewidth',1.5)
    h1 = findobj(gca,'Tag','Box');
    h2 = findobj(gca,'Tag','Lower Adjacent Value');
    h3 = findobj(gca,'Tag','Upper Adjacent Value');
    h4 = findobj(gca,'Tag','Median');
    h5 = findobj(gca,'Tag','Lower Whisker');
    h6 = findobj(gca,'Tag','Upper Whisker');
    for k=1:length(h1)
        set(h1(k),'color','k')
        set(h2(k),'color','k')
        set(h3(k),'color','k')
        set(h4(k),'color','k')
        set(h5(k),'color','k')
        set(h6(k),'color','k')
    end
    %Erosion
    sp(w(j,2)) = subplot(2,3,w(j,2));
    id3 = find(igf<-1*value);
    id4 = find(ige<-1*value);
    wbd = [igf(id3); ige(id4)];
    grp = [repmat(1,length(igf(id3)),1); repmat(2,length(ige(id4)),1)];
    plot(ones(10,1)*2.5,linspace(-0.51,0.05,10),'-k'),hold on
    bp = boxplot(wbd,grp,'symbol','.',...
        'labels',{'Flood';'Ebb'},...
        'colors',cl,...
        'boxstyle','outline',...
        'plotstyle','traditional',...
        'jitter',1,...
        'outliersize',11);
    set(bp,'linewidth',1.5)
        h1 = findobj(gca,'Tag','Box');
    h2 = findobj(gca,'Tag','Lower Adjacent Value');
    h3 = findobj(gca,'Tag','Upper Adjacent Value');
    h4 = findobj(gca,'Tag','Median');
    h5 = findobj(gca,'Tag','Lower Whisker');
    h6 = findobj(gca,'Tag','Upper Whisker');
    for k=1:length(h1)
        set(h1(k),'color','k')
        set(h2(k),'color','k')
        set(h3(k),'color','k')
        set(h4(k),'color','k')
        set(h5(k),'color','k')
        set(h6(k),'color','k')
    end
    %output to console
    disp(['Area: ' flds{j}])
    fld1 = (length([find(igf<0.001)])/(length(igf)))*100;
    fprintf('Percent flood below 1 mm: %0.2f\n',fld1)
    fld2 = (length([find(igf>-0.001)])/(length(igf)))*100;
    fprintf('Percent flood above -1 mm: %0.2f\n',fld2)
    ebb1 = (length([find(ige<0.001)])/(length(ige)))*100;
    fprintf('Percent ebb below 1 mm: %0.2f\n',ebb1)
    ebb2 = (length([find(ige>-0.001)])/(length(ige)))*100;
    fprintf('Percent ebb above -1 mm: %0.2f\n',ebb2)
end
%Plot adjustments
set([sp(1) sp(2) sp(3)],'xtick',0:0.02:0.1)
set([sp(2) sp(3) sp(5) sp(6)],'yticklabel',[])
%Top row
set(sp(1),'ylim',[-0.005 0.06],'position',[0.12 0.55 0.25 0.37])
set(sp(2),'ylim',[-0.005 0.06],'position',[0.405 0.55 0.25 0.37])
set(sp(3),'ylim',[-0.005 0.06],'position',[0.69 0.55 0.25 0.37])
%Bottom row
set(sp(4),'ylim',[-0.06 0.005],'position',[0.12 0.14 0.25 0.37])
set(sp(5),'ylim',[-0.06 0.005],'position',[0.405 0.14 0.25 0.37])
set(sp(6),'ylim',[-0.06 0.005],'position',[0.69 0.14 0.25 0.37])
title(sp(1),'Mudflat')
title(sp(2),'Fringe')
title(sp(3),'Forest')
suplabel('Net elevation change over event [mm]','y');
prettyfigures('text',12,'labels',13,'box',1)
% export_fig([sfdir 'NetChange_boxplot_v4'],'-pdf','-nocrop')
