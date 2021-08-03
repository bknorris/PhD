% Average the current meter velocity data into chunks corresponding to LL,
% ML, MH and HH tide for Paper 2 analysis
%%% Load the 2014 Data %%%%
clear
tkedir = 'd:\Projects\Mekong_F2014\DataAnalysis\Paper2\TKE\';
tkefiles = dir([tkedir '*.mat']);tkefiles = {tkefiles.name};

wvdir = 'd:\Projects\Mekong_F2014\DataAnalysis\Paper2\Spectra\';
wvfiles = dir([wvdir '*.mat']);wvfiles = {wvfiles.name};
%%%%extract TKE estimates from TKE files%%%%

%%%Extraction Settings%%%
whichwv = [1 2 4 3 4 3 5 6 7 8];
whichtke = [1 2 3 3 4 4 5 6 7 8];
ebb = [0 1 0 0 1 1 0 0 0 0];                                                    %experiments where the ebb tide was recorded
U = zeros(10,4);
V = zeros(10,4);
Uspd = zeros(10,4);
Udir = zeros(10,4);
Time = zeros(10,4);
name = cell(10,1);
inst = cell(10,1);
for i = 1:length(whichwv)
    tf = whichtke(i);load([tkedir tkefiles{tf}])
    wf = whichwv(i);load([wvdir wvfiles{wf}])
    ins = regexp(wvfiles{wf},'.+_(.*)wvstats','tokens');
    inst{i} = char(ins{:});
    name{i} = regexprep(tkefiles{tf},'TKE.mat','');
    
    %crop wvstats file to the length of the TKE record
    fn = fieldnames(Stat);
    ind = find(wvstats.time >= Stat.(fn{1}).time(1) & wvstats.time <= Stat.(fn{1}).time(end));
    
    u = wvstats.u(ind);
    v = wvstats.v(ind);
    uspd = wvstats.uspd(ind);
    udir = wvstats.udir(ind);
    wvtime = wvstats.time(ind);
    
    %%%Use times to determine the number of samples to average over by
    %dividing the timeseries into four parts (LL, ML, MH, HH)
    start = wvtime(1);stop = wvtime(end);
    [~,~,~,hr,mi,~] = datevec(stop-start);
    time = (hr*60) + mi;                                                    %calculate minutes elapsed between start/stop time of experiment
    time = datenum(0,0,0,0,floor(time/4),0);
    if ebb(i) == 0
        timee = start:time:stop;timee(end) = stop;                          %if the data is from an ebb tide, flip the time record around so the first measurement is lower tidal elevation
    elseif ebb(i) == 1
        timee = fliplr(start:time:stop);timee(1) = stop;
    end
    %%%Average into subsections of the tide:
    for j = 1:length(timee)-1
        if ebb(i) == 0
            id = find(wvtime >= timee(j) & wvtime <= timee(j+1));
        elseif ebb(i) == 1
            id = find(wvtime <= timee(j) & wvtime >= timee(j+1));
        end
        if j > 4                                                            %the FSS_1 experiment is so short that it has >4 timesteps
            continue
        else
            Uspd(i,j) = nanmean(uspd(id));
            Udir(i,j) = nanmean(udir(id));
            U(i,j) = nanmean(u(id));
            V(i,j) = nanmean(v(id));
            Time(i,j) = wvtime(id(1));
        end
    end
end
tdata.dep = name;
tdata.inst = inst;
tdata.time = Time;
tdata.u = U;
tdata.v = V;
tdata.uspd = Uspd;
tdata.udir = Udir;
clearvars -except tdata

%%%2015 Data%%%
wvdir = 'd:\Projects\Mekong_W2015\DataAnalysis\Paper2\Spectra\';
wvfiles = dir([wvdir '*.mat']);wvfiles = {wvfiles.name};
tkedir = 'd:\Projects\Mekong_W2015\DataAnalysis\Paper2\TKE\';
tkefiles = dir([tkedir '*.mat']);tkefiles = {tkefiles.name};

%wvfiles are split by experiment, not by data
whichwv = [1 2 1 2 4 5 4 5 6 7 8 8 8 9 10 10 10 10 10 10];
whichtke = [1 1 2 2 3 3 4 4 5 6 7 8 9 10 11 12 13 14 15 16];
ebb = [1 1 0 0 0 0 1 1 0 0 1 1 1 0 1 1 1 0 0 0];
U = zeros(20,4);
V = zeros(20,4);
Uspd = zeros(20,4);
Udir = zeros(20,4);
Time = zeros(20,4);
name = cell(20,1);
inst = cell(20,1);
for i = 1:length(whichwv)
    tf = whichtke(i);load([tkedir tkefiles{tf}])
    wf = whichwv(i);load([wvdir wvfiles{wf}])
    ins = regexp(wvfiles{wf},'.+_(.*)wvstats','tokens');
    inst{i} = char(ins{:});
    name{i} = regexprep(tkefiles{tf},'TKE.mat','');
    fn = fieldnames(Stat);
    
    %some experiments lasted the whole tidal cycle. Crop to low-high or
    %high-low tide
    if i >= 11 && i <= 13
        start = datenum(2015,03,09,04,30,00);
    else
        start = Stat.(fn{1}).time(1);
    end
    if i >= 5 && i <= 6
        stop = datenum(2015,03,11,16,00,00);
    elseif i == 10
        stop = datenum(2015,03,08,16,48,00);
    elseif i >= 15 && i <= 17
        stop = datenum(2015,03,13,22,00,00);
    elseif i >= 18 && i <= 20
        stop = datenum(2015,03,14,09,20,00);
    else
        stop = Stat.(fn{1}).time(end);
    end
    idx = find(Stat.(fn{1}).time >= start & Stat.(fn{1}).time <= stop);
    Stat.(fn{1}).time = Stat.(fn{1}).time(idx);
    
    %crop wvstats file to the length of the TKE record
    
    ind = find(wvstats.time >= Stat.(fn{1}).time(1) & wvstats.time <= Stat.(fn{1}).time(end));
    
    u = wvstats.u(ind);
    v = wvstats.v(ind);
    uspd = wvstats.uspd(ind);
    udir = wvstats.udir(ind);
    wvtime = wvstats.time(ind);
    
    %%%Use times to determine the number of samples to average over by
    %dividing the timeseries into four parts (LL, ML, MH, HH)
    start = wvtime(1);stop = wvtime(end);
    [~,~,~,hr,mi,~] = datevec(stop-start);
    time = (hr*60) + mi;                                                    %calculate minutes elapsed between start/stop time of experiment
    time = datenum(0,0,0,0,floor(time/4),0);
    if ebb(i) == 0
        timee = start:time:stop;timee(end) = stop;                          %if the data is from an ebb tide, flip the time record around so the first measurement is lower tidal elevation
    elseif ebb(i) == 1
        timee = fliplr(start:time:stop);timee(1) = stop;
    end
    %%%Average into subsections of the tide:
    for j = 1:length(timee)-1
        if ebb(i) == 0
            id = find(wvtime >= timee(j) & wvtime <= timee(j+1));
        elseif ebb(i) == 1
            id = find(wvtime <= timee(j) & wvtime >= timee(j+1));
        end
        if j > 4                                                            %the FSS_1 experiment is so short that it has >4 timesteps
            continue
        else
            Uspd(i,j) = nanmean(uspd(id));
            Udir(i,j) = nanmean(udir(id));
            U(i,j) = nanmean(u(id));
            V(i,j) = nanmean(v(id));
            Time(i,j) = wvtime(id(1));
        end
    end
end
tdata.dep(11:30) = name;
tdata.inst(11:30) = inst;
tdata.time(11:30,:) = Time;
tdata.u(11:30,:) = U;
tdata.v(11:30,:) = V;
tdata.uspd(11:30,:) = Uspd;
tdata.udir(11:30,:) = Udir;
clearvars -except tdata

qname = {'';'Q7';'Q5';'Q6';'Q5';'Q6';'Q4';'Q3a';'Q3b';'Q2';...
    'Q1_1A';'Q1_4A';'Q1_1A';'Q1_4A';'Q3A';'Q3B';'Q3A';'Q3B';...
    'Q2A';'Q2B';'Q2B';'Q2B';'Q2B';'Q2C';'Q4';'Q4';'Q4';'Q4';'Q4';'Q4'};
tdata.qname = qname;

%%%PLOT ROUTINE%%%
savedir = 'd:\Projects\Mekong_W2015\Figures\Paper2\CurrentDirections\';
datarep = 'd:\Projects\Documents\Writing\DataReports\';
load([datarep 'CurrentMeterPos.mat'])
load([datarep 'QuadPositions.mat'])
c = 1;
cm = brewermap(4,'YlGnBu');
colormap(cm)
for i = 2:length(qname)
    q = strcmp(tdata.qname{i},CMpos.Qname);
    p = strcmp(tdata.inst{i},CMpos.Inst);
    id = find(q == 1);id2 = find(p == 1);
    idx = intersect(id,id2);
    lat = repmat(unique(CMpos.Lat(idx)),1,4);
    lon = repmat(unique(CMpos.Lon(idx)),1,4);
    q = strcmp(tdata.qname{i},Qpos.Quad);
    id = find(q == 1);
    qlat = Qpos.Lat(id);qlon = Qpos.Long(id);

    u = tdata.u(i,:);
    v = tdata.v(i,:);
    
    
    %define map properties
    f1 = figure(c);
    set(f1,'PaperOrientation','portrait',...
        'position',[200 300   800   600]);
    set(gcf, 'PaperPositionMode', 'auto','Color','w');
    grid on
    % lon = [106.189660 106.371122];
    % lat = [9.461793 9.594033];
    arcdb = 'd:\Projects\Mekong_F2014\ArcGis\Shapefiles\';
    S = m_shaperead([arcdb 'FullCoast']);
    x = S.ncst{1}(:,1);
    y = S.ncst{1}(:,2);
    basevalue = min(y);
    p = area(x,y,basevalue);
    set(p,'FaceColor',[192,192,192]./255,...
        'EdgeColor',[250,240,230]./255,...
        'LineWidth',5);hold on
    plot(qlon,qlat,'o',...
        'Color','r',...
        'MarkerSize',18,...
        'LineWidth',1.5)    
    plot(qlon,qlat,'^',...
        'Color','k',...
        'MarkerFaceColor','k',...
        'MarkerSize',15)

    qm = zeros(4,2);
    vsc = 0.01; %vector scaling
    refscl = 0.025;
    for j = 1:4
        hold on
        qm(j) = quiver(lon(j),lat(j),vsc*u(j),vsc*v(j),0);
        set(qm(j),'Color',cm(j,:),'LineWidth',3)
    end
    qp(1) = lon(j)-0.0008;qp(2)= lat(j)-0.0009;
    qm(5) = quiver(qp(1),qp(2),vsc,0,refscl);
    text(qp(1)+5E-5,qp(2)+5E-5,[sprintf('%0.2f',refscl) ' m/s'])
    set(qm(5),'Color','k','LineWidth',3)
    hold off

    %legend
    leg = legend(qm(:,1),{'LL';'ML';'MH';'HH'});
    set(leg,'box','on',...
        'LineWidth',1.5,...
        'Color',[0.9 0.9 0.9],...
        'Position',[0.77 0.75 0.12 0.15])
    grid on
    %specify bounding box
    bbx = [lon(1)-0.001 lon(1)+0.001];
    bby = [lat(1)-0.001 lat(1)+0.001];
    set(gca,'layer','top','xlim',bbx,'ylim',bby,...
        'LineWidth',1.5,'GridLineStyle',':',...
        'FontSize',11,'FontName','Arial',...
        'TickDir','out',...
        'XMinorTick','on',...
        'YMinorTick','on',...
        'box','on',...
        'Color',[240,248,255]./255)
    xl = get(gca,'xtick');xlo = min(xl);xhi = max(xl);
    xs = xlo:0.0005:xhi;
    dms = degrees2dms(xs');
    for k = 1:length(dms)
        FRMTx{k} = [sprintf('%0.0f',dms(k,1)) sprintf('%c',char(176)) sprintf('%0.0f',dms(k,2)) sprintf('%c',char(39)) sprintf('%0.2f',dms(k,3)) '"'];
    end
    clear dms
    yl = get(gca,'ytick');ylo = min(yl);yhi = max(yl);
    ys = ylo:0.0005:yhi;
    dms = degrees2dms(ys');
    for k = 1:length(dms)
        FRMTy{k} = [sprintf('%0.0f',dms(k,1)) sprintf('%c',char(176)) sprintf('%0.0f',dms(k,2)) sprintf('%c',char(39)) sprintf('%0.2f',dms(k,3)) '"'];
    end
    set(gca,'XTick',xs,...
        'YTick',ys,...
        'XTickLabel',FRMTx,...
        'YTickLabel',FRMTy)
    title([tdata.dep{i} ' ' tdata.qname{i} ' ' tdata.inst{i}],...
        'FontSize',13,'FontName','Arial')
    xlabel('Longitude','FontSize',13,'FontName','Arial')
    ylabel('Latitude','FontSize',13,'FontName','Arial')
    c = c+1;
    fname = [tdata.dep{i} '_' tdata.qname{i} '_' tdata.inst{i}];
%     export_fig([savedir fname],'-jpeg','-nocrop')
end



