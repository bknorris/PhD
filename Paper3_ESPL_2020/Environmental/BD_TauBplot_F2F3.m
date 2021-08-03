%plot water depth, bed level, bed shear stress
clear, close all
dpath = 'e:\Mekong_W2015\DataAnalysis\Paper3\';
%First Load the Flood Experiment
tstart1 = datenum(2015,03,11,15,15,00);
tstop1 = datenum(2015,03,11,16,50,00);
bd1 = load([dpath 'BottomTrack\11-03-15\F2F3_1_bdtrace.mat']);
fn = fieldnames(bd1);
bfile1 = dir([dpath '\BedStress\11-03-15\','*.mat']);
bfile1 = {bfile1.name};
for i = 1:3
    load([dpath '\BedStress\11-03-15\' bfile1{i}]);
    bstress1.(fn{i}) = dat;
    clear dat
end
wfile1 = dir([dpath '\WaveStats\11-03-15\','*.mat']);
wfile1 = {wfile1.name};
for i = 1:3
    load([dpath '\WaveStats\11-03-15\' wfile1{i}]);
    wave1.(fn{i}) = wave;
    clear wave
end
%Then load the Ebb Experiment
tstart2 = datenum(2015,03,12,07,05,00);
tstop2 = datenum(2015,03,12,10,05,00);
bd2 = load([dpath 'BottomTrack\12-03-15\F2F3_2_bdtrace.mat']);
bfile2 = dir([dpath '\BedStress\12-03-15\','*.mat']);
bfile2 = {bfile2.name};
for i = 1:3
    load([dpath '\BedStress\12-03-15\' bfile2{i}]);
    bstress2.(fn{i}) = dat;
    clear dat
end
wfile2 = dir([dpath '\WaveStats\12-03-15\','*.mat']);
wfile2 = {wfile2.name};
for i = 1:3
    load([dpath '\WaveStats\12-03-15\' wfile2{i}]);
    wave2.(fn{i}) = wave;
    clear wave
end
%Ancillary information
d50 = [5.82E-5 1.23E-4 4.06E-5];
psand = [0.35 0.9 0.25];
pmud = [0.65 0.1 0.75];
vph = [0.061 0.061 0.061];
taucrit = [0.99 0.8724 0.299];
%Initialize Figure
f1 = figure(1);
set(f1,'PaperOrientation','portrait',...
    'position',[400 0  1000   800],...
    'renderer','painters');hold on
cl = [0.7 0.7 0.7;0.4 0.4 0.4;0 0 0];
symb = {'d','o','^'};
sp1 = zeros(3,1);
sp2 = zeros(3,1);
pl = [1 3 5;2 4 6];
for i = 1:length(fn)
    %Flood Tide
    %Plot water depth
    sp1(1) = subplot(3,2,pl(1,1));
    wvt = wave1.(fn{i}).time2;
    id = find(wvt>=tstart1&wvt<=tstop1);
    wvt = wvt(id);
    h = wave1.(fn{i}).h(id);
    plot(wvt,h,...
        'color',cl(i,:),...
        'linewidth',1.5)
    hold on
    %Plot TKE dissipation
    sp1(2) = subplot(3,2,pl(1,2));
    bdt = bd1.(fn{i}).time;
    id = find(bdt>=tstart1&bdt<=tstop1);
    bdt = bdt(id);
    bd = bd1.(fn{i}).bdist(id)-vph(i);
    plot(bdt,bd,...
        'color',cl(i,:),...
        'linewidth',1.5)
    hold on
    %Plot Bed Shear Stress
    sp1(3) = subplot(3,2,pl(1,3));
    ttime = bstress1.(fn{i}).time;
    id = find(ttime>=tstart1&ttime<=tstop1);
    ttime = ttime(id);
    tmax = bstress1.(fn{i}).tauw(id);
    plot(ttime,tmax,...
        'color',cl(i,:),...
        'linewidth',1.5)
    hold on
    taucrit = critstresswu(psand(i),pmud(i),d50(i),'cohesive');
    plot(bdt,ones(length(bdt),1)*taucrit,'-',...
        'linewidth',1.5,...
        'color',cl(i,:))
    %%%
    %Ebb Tide
    %Plot water depth
    sp2(1) = subplot(3,2,pl(2,1));
    wvt = wave2.(fn{i}).time2;
    id = find(wvt>=tstart2&wvt<=tstop2);
    wvt = wvt(id);
    h = wave2.(fn{i}).h(id);
    plot(wvt,h,...
        'color',cl(i,:),...
        'linewidth',1.5)
    hold on
    %Plot TKE dissipation
    sp2(2) = subplot(3,2,pl(2,2));
    bdt = bd2.(fn{i}).time;
    id = find(bdt>=tstart2&bdt<=tstop2);
    bdt = bdt(id);
    bd = bd2.(fn{i}).bdist(id)-vph(i);
    plot(bdt,bd,...
        'color',cl(i,:),...
        'linewidth',1.5)
    hold on
    %Plot Bed Shear Stress
    sp2(3) = subplot(3,2,pl(2,3));
    ttime = bstress2.(fn{i}).time;
    id = find(ttime>=tstart2&ttime<=tstop2);
    ttime = ttime(id);
    tmax = bstress2.(fn{i}).tauw(id);
    plot(ttime,tmax,...
        'color',cl(i,:),...
        'linewidth',1.5)
    hold on
    taucrit = critstresswu(psand(i),pmud(i),d50(i),'cohesive');
    plot(bdt,ones(length(bdt),1)*taucrit,'-',...
        'linewidth',1.5,...
        'color',cl(i,:))
end
t1 = bstress1.(fn{i}).time;
t2 = bstress2.(fn{i}).time;
set(sp1,'xlim',[t1(1) t1(end)])
set(sp2,'xlim',[t2(1) t2(end)],'yticklabel',[])
set([sp1(1) sp1(2) sp2(1) sp2(2)],'xticklabel',[])
set(sp1(1),'position',[0.12 0.7 0.4 0.22],'ylim',[0 1.5])
set(sp1(2),'position',[0.12 0.42 0.4 0.22],'ylim',[-0.04 0.04])
set(sp1(3),'position',[0.12 0.15 0.4 0.22],'ylim',[0 2.5])
set(sp2(1),'position',[0.55 0.7 0.4 0.22],'ylim',[0 1.5])
set(sp2(2),'position',[0.55 0.42 0.4 0.22],'ylim',[-0.04 0.04])
set(sp2(3),'position',[0.55 0.15 0.4 0.22],'ylim',[0 2.5])
datetick(sp1(3),'x','HH:MM','keepticks','keeplimits')
datetick(sp2(3),'x','HH:MM','keepticks','keeplimits')
xlabel(sp1(3),['Time on ' datestr(t1(1),'dd-mm-yy')])
xlabel(sp2(3),['Time on ' datestr(t2(1),'dd-mm-yy')])
title(sp1(1),'Flood Tide')
title(sp2(1),'Ebb Tide')
ylabel(sp1(1),'Water Depth [m]')
ylabel(sp1(2),'BLE [m]')
ylabel(sp1(3),'\tau_b [Pa]')
prettyfigures('text',12,'labels',13,'box',1)
sfdir = 'd:\GradSchool\DataAnalysis\Paper3\Figures\';
export_fig([sfdir 'BstressExamp_F2F3'],'-pdf','-nocrop')