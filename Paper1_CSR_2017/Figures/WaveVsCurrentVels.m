%Routine to plot a new figure for Paper 1. This figure will try to
%differentiate between current and wave-induced turbulence. Here, I will
%find wave velocities that exceed 2*mag of mean current. These velocities
%will be plotted in one way, and those of the current will be plotted in
%another. 

%make sure 'CatTurbData_v2.mat' is in the path

clear
%load the vegdat files and name based on the size and stage of tide
vpdir = 'd:\Projects\Mekong_W2015\DataAnalysis\DataReports\Paper1_QuadAnalysis\ThirdAttempt\';
vd = struct(); %big vegdat structure
order = [2 4 3 1]; %load files in order

%%%10cm%%%
files = dir([vpdir '*local10cm*.mat']);files = {files.name};
for i = 1:length(files)
    load([vpdir files{order(i)}])
    stage = regexp(files{order(i)},'.+_(.*).mat','tokens');
    fname = ['ten' char(stage{:})];
    vd.(fname).n = vegdat.n;
    vd.(fname).meanD = vegdat.meanD;
    vd.(fname).a = vegdat.a;
    vd.(fname).phi = vegdat.phi;
end
%%%20cm%%%
files = dir([vpdir '*local20cm*.mat']);files = {files.name};
for i = 1:length(files)
    load([vpdir files{order(i)}])
    stage = regexp(files{order(i)},'.+_(.*).mat','tokens');
    fname = ['twe' char(stage{:})];
    vd.(fname).n = vegdat.n;
    vd.(fname).meanD = vegdat.meanD;
    vd.(fname).a = vegdat.a;
    vd.(fname).phi = vegdat.phi;
end
%%%30cm%%%
files = dir([vpdir '*local30cm*.mat']);files = {files.name};
for i = 1:length(files)
    load([vpdir files{order(i)}])
    stage = regexp(files{order(i)},'.+_(.*).mat','tokens');
    fname = ['thr' char(stage{:})];
    vd.(fname).n = vegdat.n;
    vd.(fname).meanD = vegdat.meanD;
    vd.(fname).a = vegdat.a;
    vd.(fname).phi = vegdat.phi;
end
%%%40cm%%%
files = dir([vpdir '*local40cm*.mat']);files = {files.name};
for i = 1:length(files)
    load([vpdir files{order(i)}])
    stage = regexp(files{order(i)},'.+_(.*).mat','tokens');
    fname = ['fou' char(stage{:})];
    vd.(fname).n = vegdat.n;
    vd.(fname).meanD = vegdat.meanD;
    vd.(fname).a = vegdat.a;
    vd.(fname).phi = vegdat.phi;
end
%%%50cm%%%
files = dir([vpdir '*local50cm*.mat']);files = {files.name};
for i = 1:length(files)
    load([vpdir files{order(i)}])
    stage = regexp(files{order(i)},'.+_(.*).mat','tokens');
    fname = ['fif' char(stage{:})];
    vd.(fname).n = vegdat.n;
    vd.(fname).meanD = vegdat.meanD;
    vd.(fname).a = vegdat.a;
    vd.(fname).phi = vegdat.phi;
end
%%%1m%%%
load([vpdir 'Vegdat_1m.mat'])
vd.full.n = vegdat.n;
vd.full.meanD = vegdat.meanD;
vd.full.a = vegdat.a;
vd.full.phi = vegdat.phi;
vd.X = vegdat.Xshore;
clear vegdat
vegdat = vd;

%%%Velocity Statistics%%%
CatTurbData_v2
fn = fieldnames(vegdat);
E = [veldat.four.E; veldat.five.E];
std = [veldat.four.Estd; veldat.five.Estd];
H = [veldat.four.depth; veldat.five.depth];
Hs = [veldat.four.wrms; veldat.five.wrms].*sqrt(2);
Uc = [veldat.four.uc; veldat.five.uc];
Uw = [veldat.four.uw; veldat.five.uw];
Uwstd = [veldat.four.uwstd; veldat.five.uwstd];
Umag = [veldat.four.umag; veldat.five.umag];
gamma = Hs./H;

%To differentiate wave versus current-induced turbulence, find values of
%2*Uc that are less than the wave velocity, Uw
f1 = figure(1);
set(f1,'PaperOrientation','portrait',...
    'position',[400 100   800 600]);
c = brewermap(4,'YlGnBu');colormap(c)
c(1,:) = [255 235 0]./255;
vid = strfind(fn,'twe');
id = find(not(cellfun('isempty',vid)));
a = zeros(54,4);
for i = 1:length(id)
    a(:,i) = vegdat.(fn{id(i)}).a;
end
tide = {'LL';'ML';'MH';'HH'};
for i = 1:4
    id = find(2.*Uc(:,i) > Uw(:,i));id2 = setxor(1:54,id);
    tkew = E(id2,i);
    tkec = E(id,i);
    aw = a(id2,i);
    ac = a(id,i);
    zid = find(aw>0);aw = aw(zid);tkew = tkew(zid);
    zid = find(ac>0);
    if isempty(zid)
        ac = ac;tkec = tkec;
    else
        ac = ac(zid);tkec = tkec(zid);
    end
    p(1) = plot(aw,tkew,'^','color','k',...
        'linewidth',1.5,...
        'markerfacecolor',c(i,:),...
        'markersize',10); hold on
    p(2) = plot(ac,tkec,'o','color','k',...
        'linewidth',1.5,...
        'markerfacecolor',c(i,:),...
        'markersize',10);
    disp(tide{i})
    disp(['Percent wave-dominated: ' num2str((length(aw)/(length(aw)+length(ac)))*100) '%'])
    disp(['Percent current-dominated: ' num2str((length(ac)/(length(aw)+length(ac)))*100) '%'])
end
legend(p,{'Wave-dominated';'Current-dominated'})
title('wave vs current turbulence by vegetation density')
ylabel('\epsilon (Wkg^-^1)')
xlabel('a (m^-^1)')
set(gca,'yscale','log')
prettyfigures('box',1)


f2 = figure(2);
set(f2,'PaperOrientation','portrait',...
    'position',[400 100   1000 400]);
set(gcf,'color','w','PaperPositionMode','auto')
lins = {':','-.','--','-'};

xs = -60:10:100;
X = vegdat.X;

for i = 1:4
    Ez = E(:,i);
    id = find(2.*Uc(:,i) > Uw(:,i));id2 = setxor(1:54,id);
    Ew = Ez(id2);
    xs = X(id2);

    pp(1) = plot(xs,Ew,'^','Color','k',...
        'LineWidth',1.5,...
        'markerfacecolor',c(i,:),...
        'markeredgecolor','k',...
        'markersize',10); hold on
    
    Ec = Ez(id);
    xs = X(id);
    pp(2) = plot(xs,Ec,'o','Color','k',...
        'LineWidth',1.5,...
        'markerfacecolor',c(i,:),...
        'markeredgecolor','k',...
        'markersize',10); hold on
end
legend(pp,{'Wave-dominated';'Current-dominated'})
set(gca,'LineWidth',1.5,...
    'FontSize',14,...
    'FontName','Arial',...
    'TickDir','out',...
    'YScale','log')
ylabel('\epsilon (Wkg^-^1)','FontSize',18)
xlabel('Cross-shore distance (m)','FontSize',18)
title('Wave vs Current Turbulence, by x-shore distance')

