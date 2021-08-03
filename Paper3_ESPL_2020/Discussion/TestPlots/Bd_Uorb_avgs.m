clear, close all
load('g:\Mekong_W2015\DataAnalysis\Paper3\Wavelet\11-03-15\waveIGevents.mat')
sfdir = 'f:\GradSchool\DataAnalysis\Paper3\WorkingFigures\Wavelet\';
expname = 'F2F3_1';
exploc = {'Mudflat';'Fringe';'Forest'};
% exploc = {'x = -10cm';'x = 10cm';'x = 20cm'};

ltxt = cell(3,1);
p1 = zeros(3,1);
fn = fieldnames(eventd);
cl = [207 176 126;
    60 166 74;
    4 76 41]./255;
symb = {'o';'^';'s'};
f1 = figure(1);
set(f1,'PaperOrientation','portrait',...
    'position',[400 100   800   400],...
    'renderer','painters')
sp = zeros(3,1);
% spn = [1 2;3 4;5 6];
for i = 1:length(fn);
%     sp(i) = subplot(3,1,i);
    period = [eventd.(fn{i}).wave.period eventd.(fn{i}).ig.period];
    bdvar = [eventd.(fn{i}).wave.bdvar; eventd.(fn{i}).ig.bdvar];
    uorb = [eventd.(fn{i}).wave.orbwv; eventd.(fn{i}).ig.orbwv];
    uorb = uorb(~isnan(uorb));
    bdvar = bdvar(~isnan(bdvar));
    uorb = uorb(1:length(bdvar));
    bins = linspace(min(uorb),max(uorb),30);
    %%%
    [b,n,s] = bindata(uorb,bdvar,bins);
    p1(i) = errorbar(bins,b,s,symb{i},...
        'color',cl(i,:),...
        'linewidth',1.5,...
        'markersize',8,...
        'markeredgecolor','k',...
        'markerfacecolor',cl(i,:));hold on

end
% set(gca,'yscale','log')
% leg = legend(p1(p1~=0),ltxt(~cellfun('isempty',ltxt)));
%Global plot adjustments
break
set(sp(1),'ylim',[0 25E-6],'xlim',[0.1 0.5],...
    'position',[0.12 0.71 0.75 0.22],...
    'xticklabel',[])
set(sp(2),'ylim',[0 4E-8],'xlim',[0.1 0.5],...
    'position',[0.12 0.41 0.75 0.22],...
    'xticklabel',[])
set(sp(3),'ylim',[0 5E-7],'xlim',[0.1 0.5],...
    'position',[0.12 0.11 0.75 0.22])
ylabel(sp(1),'\sigma_{bd}/\Deltat*U [m]')
ylabel(sp(2),'\sigma_{bd}/\Deltat*U [m]')
ylabel(sp(3),'\sigma_{bd}/\Deltat*U [m]')
xlabel(sp(3),'Representative wave orbital velocity [ms^{-1}]')
% prettyfigures('text',12,'labels',13,'box',1)
% set(leg,'position',[0.61,0.17,0.1 0.2])
% export_fig([sfdir expname '_evlBDavg_acc_ero'],'-png')
%
