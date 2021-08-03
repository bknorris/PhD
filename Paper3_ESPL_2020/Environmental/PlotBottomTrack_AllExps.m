%Plot all bed level traces for Paper 3
clear,close all
dpath = '\DataAnalysis\Paper3\';
ypath = 'Mekong_W2015';
sdir = 'd:\Mekong_W2015\DataAnalysis\Paper3\Turbulence\Figures\';
%load run file to tell program which files & VPs to load
rdir = 'd:\MemoryStick\GradSchool\DataAnalysis\Paper3\ExperimentalDesign\';
fid = fopen([rdir 'AutoBedLevelPlot_v2.csv']);
rfile = textscan(fid,'%s%s%s%n%s','delimiter',',');
rdate = rfile{1}; %date of expt.
rarea = rfile{2}; %which area [mud/fringe/forest]
rtide = rfile{3}; %which tide [flood/ebb]
rvph = rfile{4}; %vecpro height above bed
rht = rfile{5}; %approximate high tide time
cc = jet(length(rvph));
npath = ['d:\' ypath dpath];
folders = dir([npath '\VPs\']);
folders = {folders(3:end).name};
%% Initialize Figure
f1 = figure(1);
set(f1,'PaperOrientation','portrait',...
    'position',[400 100   800   500]);
sp(1) = subplot(311);hold on
plot(zeros(1,10),linspace(-40,40,10),'-k','linewidth',1.5);
sp(2) = subplot(312);hold on
plot(zeros(1,10),linspace(-40,40,10),'-k','linewidth',1.5);
sp(3) = subplot(313);hold on
plot(zeros(1,10),linspace(-40,40,10),'-k','linewidth',1.5);
%% Plot bd traces
for ii = 1:length(folders)
    if any(strcmp(folders{ii},rdate))
        vpid = strcmp(folders{ii},rdate);vpid = find(vpid);
        disp(['Loading files from ' npath 'BottomTrack\' folders{ii} '\'])
        fname = dir([npath 'BottomTrack\' folders{ii} '\' '*bdtrace.mat']);
        bd = load([npath 'BottomTrack\' folders{ii} '\',fname.name]);
    else
        continue
    end
    for j = 1:3
        dfn = fieldnames(bd);
        bdstart = find(~isnan(bd.(dfn{j}).bdist),1,'first');
        bdlvl = (bd.(dfn{j}).bdist(bdstart)-bd.(dfn{j}).bdist)*1000;
        bdt = bd.(dfn{j}).time;
        nid = find(isnan(bdlvl));
        bdt(nid) = [];
        bdlvl(nid) = [];
%         bdlvl = fastsmooth(bdlvl,400,1,1);
        %Find time till/from high tide
        hightide = datenum([rdate{vpid(j)} ' ' rht{vpid(j)}],'dd-mm-yy HH:MM:SS');
        tn = length(bdt);
        toffset(1) = etime(datevec(bdt(1)),datevec(hightide))/60; %time in minutes
        toffset(2) = etime(datevec(bdt(end)),datevec(hightide))/60;
        timeadj = linspace(toffset(1),toffset(2),tn);
        if strcmp(rtide(vpid(j)),'flood')
            timeadj(timeadj>0) = NaN; %remove values after HT
        elseif strcmp(rtide(vpid(j)),'ebb')
            timeadj(timeadj<0) = NaN; %remove values before HT
        end
        if strcmp(rarea(vpid(j)),'mud')
            sp(1) = subplot(311);
            plot(timeadj,bdlvl,'color',cc(vpid(j),:),'linewidth',1.5)
        elseif strcmp(rarea(vpid(j)),'fringe')
            sp(2) = subplot(312);
            plot(timeadj,bdlvl,'color',cc(vpid(j),:),'linewidth',1.5)
        else
            sp(3) = subplot(313);
            plot(timeadj,bdlvl,'color',cc(vpid(j),:),'linewidth',1.5)
        end
    end
end
%stick in the bdtrace from 09-03-15 bc it has a different data format
load('e:\Mekong_W2015\DataAnalysis\Paper3\BottomTrack\BDtrace_floodebb.mat')
hightide = datenum(2015,03,09,04,40,00);
bdt = data.vp3.time;bdlvl = data.vp3.ble*1000;
toffset(1) = etime(datevec(bdt(1)),datevec(hightide))/60; %time in minutes
toffset(2) = etime(datevec(bdt(end)),datevec(hightide))/60;
tn = length(bdt);
timeadj = linspace(toffset(1),toffset(2),tn);
flood = find(timeadj<0);ebb = find(timeadj>0);
sp(2) = subplot(312);
plot(timeadj(flood),bdlvl(flood),'m','linewidth',1.5)
plot(timeadj(ebb),bdlvl(ebb),'k','linewidth',1.5)

%% Global plot adjustments
set(sp,'xlim',[-150 175],'ylim',[-40 40])
xlabel(sp(3),'Minutes to/from High Tide')
ylabel(sp(1),'Bed Level [mm]')
ylabel(sp(2),'Bed Level [mm]')
ylabel(sp(3),'Bed Level [mm]')
prettyfigures('text',13,'labels',14,'box',1)
% export_fig([sdir 'Btrace_allExps_v2'],'-pdf','-nocrop')