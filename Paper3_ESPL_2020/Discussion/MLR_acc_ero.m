%Testbed
clear, close all
load('d:\Mekong_W2015\DataAnalysis\Paper3\Wavelet\CombEventData_flood.mat');
flds = fieldnames(dat);
vars = {'depth','usqd','taub','eps';...
    'depth','orbwv','usqd','eps';...
    'depth','taub','eps',NaN};
beta = [0.25 -0.04 -0.62 0.14;...
    0.26 -0.46 -0.17 -0.20;...
    0.23 -0.05 -0.48 0.51];
dbd = [6.68;-0.48;9.95];
symb = {'s';'o';'d'};
cl = [0.6 0.6 0.6;0.4 0.4 0.4;0 0 0];
figure(), hold on
for i = 1:length(flds)
    bdmed = [dat.(flds{i}).wave.bdmed; dat.(flds{i}).ig.bdmed];
    a = [dat.(flds{i}).wave.(vars{i,1}); dat.(flds{i}).ig.(vars{i,1})];
    b = [dat.(flds{i}).wave.(vars{i,2}); dat.(flds{i}).ig.(vars{i,2})];
    c = [dat.(flds{i}).wave.(vars{i,3}); dat.(flds{i}).ig.(vars{i,3})];
    if ischar(vars{i,4})
        d = [dat.(flds{i}).wave.(vars{i,4}); dat.(flds{i}).ig.(vars{i,4})];
    else
        d = 0;
    end
    xs = beta(i,1)*max(a)+beta(i,2)*max(b)+beta(i,3)*max(c)+beta(i,4)*max(d);
    y = dbd(i);
    plot(abs(xs),y,symb{i},...
        'markerfacecolor',cl(i,:),...
        'markeredgecolor','k',...
        'linewidth',1.5,...
        'markersize',8)

end