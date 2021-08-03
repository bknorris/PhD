%Plot bed shear stress by the net bed change
clear, close all
bd = load('G:\Mekong_W2015\DataAnalysis\Paper3\BottomTrack\05-03-15\F2F2_1_bdtrace.mat');
load('D:\Projects\Mekong_W2015\DataAnalysis\Paper3\RS\Two\F2F2_3_RS10min.mat');
load('D:\Projects\Mekong_W2015\DataAnalysis\Paper1\AveragedVelocities\Two\F2F2_1Vave.mat')
fn = fieldnames(bd);
f1 = figure(1);
set(f1,'PaperOrientation','portrait',...
    'position',[400 100   500   400],...
    'renderer','painters');hold on
cl = [207 176 126;
    60 166 74;
    4 76 41]./255;
rho = 1029; %kg/m^3
p = zeros(3,1);
bins = [25 25 1];
for i = 1:3
    %bed shear stress
    uw = RS.f2f2.(fn{i}).uw./100;
    taub = abs(uw(:,bins(i))*rho);
    %bottom variance
    t1 = Avgs.(fn{i}).time;
    t2 = bd.(fn{i}).time;
    bdvar = zeros(length(t1),1);
    for ii = 1:length(t1)-1
        id = find(t2>=t1(ii)&t2<=t1(ii+1));
        bdvar(ii) = var(bd.(fn{i}).bdist(id));
    end
    U = nanmean(Avgs.(fn{i}).Umag(1:25,:));
    bdn = bdvar./(U'*300);
%     bdn(isnan(bdn)) = 0;
%     bdn(isnan(taub)) = [];
%     taub(isnan(taub)) = [];
%     taub(bdn==0) = [];
%     bdn(bdn==0) = [];
%     pf = polyfit(log10(taub),log10(bdn),1);
%     xs = linspace(min(taub),max(taub),length(taub));
%     y_hat = exp(pf(1)*log(xs)+pf(2));
%     plot(xs,y_hat,'color',cl(i,:),...
%         'linewidth',1.5)
    p(i) = plot(taub,bdn,'*',...
        'color',cl(i,:),...
        'markersize',6); hold on
end
set(gca,'yscale','log')