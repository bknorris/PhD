%Plot var(bd) versus epsilon for selected experiments. Look for positive
%correlation.
clear, close all
bd = load('g:\Mekong_W2015\DataAnalysis\Paper3\BottomTrack\05-03-15\F2F2_1_bdtrace.mat');
load('g:\Mekong_W2015\DataAnalysis\Paper3\Turbulence\V1\5March2015_VelsTKE.mat')
load('g:\Mekong_W2015\DataAnalysis\Paper1\AveragedVelocities\Two\F2F2_1Vave.mat')
fn = fieldnames(bd);
f1 = figure(1);
set(f1,'PaperOrientation','portrait',...
    'position',[400 100   500   400],...
    'renderer','painters');hold on
symb = {'o';'^';'s'};
cl = [0.7 0.7 0.7;0.4 0.4 0.4;0 0 0];
a = [0,1.86,0.29];
for i = 1:length(fn);
    t1 = Stat.(fn{i}).time;
    t2 = bd.(fn{i}).time;
    bdvar = zeros(length(t1),1);
    for ii = 1:length(t1)-1
        id = find(t2>=t1(ii)&t2<=t1(ii+1));
        bdvar(ii) = var(bd.(fn{i}).bdist(id));
    end
    U = nanmean(Avgs.(fn{i}).Umag(1:25,:));
    nid = find(~isnan(bd.(fn{i}).bdist),1,'first');
    bdn = bdvar./(U'*300);
%     eps = nanmean((Stat.(fn{i}).z1.E+Stat.(fn{i}).z2.E)./2);
%     plot(eps(1:end-1),bdn(1:end-1),symb{i},...
%         'markerfacecolor',cl(i,:),...
%         'markeredgecolor','k',...
%         'markersize',8)
      plot(a(i),nanmean(bdn(1:end-1)),symb{i},...
        'markerfacecolor',cl(i,:),...
        'markeredgecolor','k',...
        'markersize',8)
end
% set(gca,'yscale','log','xscale','log')