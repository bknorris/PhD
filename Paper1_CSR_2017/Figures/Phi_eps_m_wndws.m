%Plot routine to generate cross-shore, phi/number of stems/stem diameter
%plots, plus phi,epsilon plots for various tide stages.

%make sure 'CatTurbData_v2.mat' is in the path

clear
%load the vegdat files and name based on the size and stage of tide
vpdir = 'd:\Projects\Documents\Writing\DataReports\Volumetric\';
vd = struct(); %big vegdat structure
order = [2 4 3 1]; %load files in order

%%%10cm%%%
files = dir([vpdir '*local10cm*.mat']);files = {files.name};
for i = 1:length(files)
    load([vpdir files{order(i)}])
    stage = regexp(files{order(i)},'.+_(.*)_vol.mat','tokens');
    fname = ['ten' char(stage{:})];
    vd.(fname).n = vegdat.n;
    vd.(fname).meanD = vegdat.meanD;
    vd.(fname).a = vegdat.a_vol;
    vd.(fname).phi = vegdat.phi_vol;
end
%%%20cm%%%
files = dir([vpdir '*local20cm*.mat']);files = {files.name};
for i = 1:length(files)
    load([vpdir files{order(i)}])
    stage = regexp(files{order(i)},'.+_(.*)_vol.mat','tokens');
    fname = ['twe' char(stage{:})];
    vd.(fname).n = vegdat.n;
    vd.(fname).meanD = vegdat.meanD;
    vd.(fname).a = vegdat.a_vol;
    vd.(fname).phi = vegdat.phi_vol;
end
%%%30cm%%%
files = dir([vpdir '*local30cm*.mat']);files = {files.name};
for i = 1:length(files)
    load([vpdir files{order(i)}])
    stage = regexp(files{order(i)},'.+_(.*)_vol.mat','tokens');
    fname = ['thr' char(stage{:})];
    vd.(fname).n = vegdat.n;
    vd.(fname).meanD = vegdat.meanD;
    vd.(fname).a = vegdat.a_vol;
    vd.(fname).phi = vegdat.phi_vol;
end
%%%40cm%%%
files = dir([vpdir '*local40cm*.mat']);files = {files.name};
for i = 1:length(files)
    load([vpdir files{order(i)}])
    stage = regexp(files{order(i)},'.+_(.*)_vol.mat','tokens');
    fname = ['fou' char(stage{:})];
    vd.(fname).n = vegdat.n;
    vd.(fname).meanD = vegdat.meanD;
    vd.(fname).a = vegdat.a_vol;
    vd.(fname).phi = vegdat.phi_vol;
end
%%%50cm%%%
files = dir([vpdir '*local50cm*.mat']);files = {files.name};
for i = 1:length(files)
    load([vpdir files{order(i)}])
    stage = regexp(files{order(i)},'.+_(.*)_vol.mat','tokens');
    fname = ['fif' char(stage{:})];
    vd.(fname).n = vegdat.n;
    vd.(fname).meanD = vegdat.meanD;
    vd.(fname).a = vegdat.a_vol;
    vd.(fname).phi = vegdat.phi_vol;
end
%%%1m%%%
load([vpdir 'Vegdat_local1m_vol.mat'])
vd.full.n = vegdat.n;
vd.full.meanD = vegdat.meanD;
vd.full.a = vegdat.a_vol;
vd.full.phi = vegdat.phi_vol;
vd.X = vegdat.Xshore;
clear vegdat
vegdat = vd;

%%%Velocity Statistics%%%
CatTurbData_v2
fn = fieldnames(vegdat);
E = [veldat.four.E; veldat.five.E];
stdv = [veldat.four.Estd; veldat.five.Estd];
H = [veldat.four.depth; veldat.five.depth];
Hs = [veldat.four.wrms; veldat.five.wrms].*sqrt(2);
gamma = Hs./H;
savefigdir = 'd:\Projects\Mekong_W2015\Figures\Paper2\Phi-Eps\';


%%%Plot Routine%%%
f1 = figure(1);
set(f1,'PaperOrientation','portrait',...
    'position',[400 100   1000 800]);
set(gcf,'color','w','PaperPositionMode','auto')
tt = {'LL','ML','MH','HH'};
cb = {'YlOrRd';'PuBu';'BuPu';'YlGn'};
symb = {'^';'d';'s';'o';'p'};
sp = zeros(1,4);
p = zeros(1,5);

for i = 1:4
    c = brewermap(5,cb{i});
    vid = strfind(fn,tt{i});
    id = find(not(cellfun('isempty',vid)));
    sp(i) = subplot(2,2,i);
    for j = 1:length(id)
            phi = vegdat.(fn{id(j)}).phi;
            zid = find(phi == 0);
            zid = setxor(1:54,zid);
            Ez = E(zid,i);Ez(Ez > 1E-3) = NaN;
            phiz = phi(zid);
            bins = linspace(min(phiz),max(phiz),5);

            [b,n,s] = bindata(phiz,Ez,bins);
        plot(bins,b,'-','LineWidth',1.5,...
            'Color',c(j,:)),hold on
%         plot(phimu,epsmu,'o','LineWidth',1.5,...
%             'Color',c(j,:)),hold on
        p(j) = errorbar(bins,b,s,symb{j},...
            'Color','k',...
            'LineWidth',1.5,...
            'markerfacecolor',c(j,:),...
            'markeredgecolor','k',...
            'markersize',10);
    end
end
