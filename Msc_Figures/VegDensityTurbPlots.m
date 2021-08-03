%Plot routine to generate figures for the short-format paper.

%make sure 'CatTurbData_v2.mat' is in the path

clear
%load the vegdat files and name based on the size and stage of tide
vpdir = 'd:\Projects\Documents\Writing\DataReports\ThirdAttempt\';
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
ucubed = [veldat.four.ucubed; veldat.five.ucubed];
gamma = Hs./H;
savefigdir = 'd:\Projects\Mekong_W2015\Figures\Paper2\Phi-Eps\';

%%%Plot Routine%%%
symb = {'^','>','d','s','<'};
tt = {'ten';'twe';'thr';'fou';'fif'};
c = [101 45 144;182 36 103;...
    240 90 34;252 184 19;109 199 222]./255;

f1 = figure(1);
set(f1,'PaperOrientation','portrait',...
    'position',[400 100   1000 800]);
set(gcf,'color','w','PaperPositionMode','auto')
%%%%
for i = 1:5
sp(i) = subplot(1,5,i);
vid = strfind(fn,tt{i});
id = find(not(cellfun('isempty',vid)));
sphi = zeros(54,4);sa = zeros(54,4); %'s' for 'scale'
for j = 1:length(id)
    sphi(:,j) = vegdat.(fn{id(j)}).phi;
    sa(:,j) = vegdat.(fn{id(j)}).a;
end
idx = find(sa == 0);
sa(idx) = NaN;
fa = vegdat.full.a;
emean = mean(E,2);
samean = mean(sa.*ucubed,2);
plot(fa.*(mean(ucubed,2)),emean,'o',...
    'MarkerSize',5,...
    'MarkerFaceColor',[0.8 0.8 0.8],...
    'MarkerEdgeColor','k'), hold on
plot(samean,emean,symb{i},...
    'MarkerSize',10,...
    'MarkerFaceColor',c(i,:),...
    'MarkerEdgeColor','k')
%%%%
%linreg
samean(isnan(samean)) = 0;
pf = polyfit(samean,emean,1);
stats = regstats(emean,samean,'linear','fstat');
disp([tt{i} ' line slope ' num2str(pf(1))])
disp([tt{i} ' p-value ' num2str(stats.fstat.pval)])
pv = polyval(pf,samean);
resid = emean-pv;
SSresid = sum(resid.^2);
SStotal = (length(emean)-1)*var(emean);
rsq = 1-SSresid/SStotal;
disp([tt{i} ' rsq ' num2str(rsq)])

set(gca,'LineWidth',1.5,...
    'FontSize',14,...
    'FontName','Arial',...
    'XScale','log',...
    'YScale','log',...
    'TickDir','out',...
    'FontSize',14,...
    'FontName','Arial')
end
set(sp(2),'position',[0.55 0.1 0.4 0.85])
set(sp(1),'position',[0.1 0.76 0.3 0.19],'XTickLabel',[])
set(sp(3),'position',[0.1 0.54 0.3 0.19],'XTickLabel',[])
set(sp(4),'position',[0.1 0.32 0.3 0.19],'XTickLabel',[])
set(sp(5),'position',[0.1 0.1 0.3 0.19])
set(sp,'XLim',[1E-6 1E-1],'YLim',[1E-6 1E-2])
ylabel(sp(2),'\epsilon (Wkg^-^1)')
xlabel(sp(2),'a|u|^3')
ylabel(sp(3),'\epsilon (Wkg^-^1)')
xlabel(sp(5),'a|u|^3')

pf = polyfit(fa.*(mean(ucubed,2)),emean,1);
stats = regstats(emean,fa.*(mean(ucubed,2)),'linear','fstat');
disp(['full line slope ' num2str(pf(1))])
disp(['full p-value ' num2str(stats.fstat.pval)])
pv = polyval(pf,fa.*(mean(ucubed,2)));
resid = emean-pv;
SSresid = sum(resid.^2);
SStotal = (length(emean)-1)*var(emean);
rsq = 1-SSresid/SStotal;
disp(['full rsq ' num2str(rsq)])