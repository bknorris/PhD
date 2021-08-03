%Plot bed aggredation or scour for a given HTA or VTA deployment

clear
datdir = 'D:\Projects\Mekong_W2015\Data\Vectrino\7March2015\Vectrinos\';
datdir2 = 'D:\Projects\Mekong_W2015\DataAnalysis\VPs\Paper1\';
figdir = 'D:\Projects\Mekong_W2015\Figures\Environment\';

fname = dir([datdir '*070315.mat']);
fname = {fname(1:3).name};
hab = [0.063 0.061 0.062]*100;
% hab = [0.242 0.271 0.24]*100;
%load the data. this step takes the longest.
name = {'v1';'v2';'v3'};
for i = 1:3
    load([datdir fname{i}])
    dat.(name{i}).bdist = VPRO.Data.BottomCheck_BottomDistance;
    bchtm = VPRO.Data.BottomCheck_HostTimeMatlab;
    gmt2ict = datenum(0,0,0,1,0,0)*7;
    dat.(name{i}).bchtm = bchtm+gmt2ict;
    clear VPRO
end
save([datdir2 'HTAday1BdTrack'],'dat','-v7.3')
f = figure;
set(f,'PaperOrientation','portrait',...
    'position',[400 200   600   500]);
set(gcf,'color','w','PaperPositionMode','auto')
markers = {'^','s','o'};
c = [0 0 0;1 0 0;0 0 1];
for i = 1:3
    h = hab(i);
    bdist = dat.(name{i}).bdist.*100;
    bchtm = dat.(name{i}).bchtm;
    bdmax = runningmax(bdist,100);
    bdmax = my_running_median(bdmax,500); %despike
    bdmax = smooth(bdmax,100,'sgolay');

    %window to smooth the line further
    win = 400; 
    step = 200;
    ind = [1 step:step:length(bdmax)];
    bdavs = zeros(length(ind),1);
    bdavt = zeros(length(ind),1);
    for ii = 1:length(ind)
        if abs(length(bdmax)-ind(ii)) < win
            continue
        else
            bwin = bdmax(ind(ii):ind(ii)+win);
            bdavs(ii,:) = nanmean(bwin);
            bdavt(ii,:) = bchtm(ind(ii));
        end
    end
    bdavs(bdavs == 0) = []; %remove trailing zeros
    bdavt(bdavt == 0) = [];
    sed = h-bdavs;
    p(i) = line_fewer_markers(bdavt,sed,15,['-' markers{i}],'Color',c(i,:),'MarkerSize',8,...
        'LineWidth',1.5);hold on
    
end
datetick('x','HH:MM','keepticks','keeplimits')

%draw refline about zero
rlx = linspace(bdavt(1),bdavt(end),length(bdavt));rly = zeros(1,length(rlx));
rl = plot(rlx,rly,'--k','LineWidth',1.5);hold off

%legend
leg = legend(p,{'HTA VP1';'HTA VP2';'HTA VP3'});
set(leg,'box','off','position',[0.25 0.8 0.05 0.05],...
'FontAngle','italic')%'FontWeight','Bold')

%global adjustments
set(gca,'Box','on','LineWidth',1.5,'FontSize',12,...
    'Xlim',[bdavt(1) bdavt(end)],'Ylim',[-2 2])
xlabel('\itTime on 07/03/2015')
ylabel('\itBed Level Change (cm)')

% title('Bottom Tracking, HTA Day 1')
% export_fig([figdir 'HTAday1BTrack'],'-jpeg')

