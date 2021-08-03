%Plot TKE and BLE time series for the experiment on 09-03-15
%Then plot TKE vs. BLE for the same experiment
clear
load('E:\Mekong_W2015\DataAnalysis\Paper3\BottomTrack\BDtrace_floodebb.mat')
tke = load('E:\Mekong_W2015\DataAnalysis\Paper3\Turbulence\09-03-15\HTA_09_TKE_v3.mat');
tfn = fieldnames(tke);
bfn = fieldnames(data);
f(1) = figure(1);
c = [207 176 126;
    60 166 74;
    4 76 41]./255;
for j = 1:length(tfn)
    bh = data.(bfn{j}).ble;
    bdt = data.(bfn{j}).time;
    tdt = tke.(tfn{j}).time;
    eps = (tke.(tfn{j}).z1.E+tke.(tfn{j}).z2.E)./2;
    [~,idx] = unique(tdt);
    ep10 = zeros(30,length(bdt));
    for k = 1:30
        ep10(k,:) = interp1(tdt(idx),eps(k,idx),bdt,'linear','extrap');
    end
    %04/17/20: Find depth bins that are a consistent elevation
    %above bed level for TKE extraction!
    bstart = find(bh,1,'first');
    vh = 0.24-0.04-linspace(0.001,0.03,30);
    eps = zeros(length(bdt),1);bb = zeros(length(bdt),1);
    for k = 1:length(bdt)
%         [minval,bid] = min(abs((bh(k)+0.007)-vh));
        [~,bid] = min(abs((bh(k)+0.04)-vh));
        eps(k) = ep10(bid,k);
        bb(k) = vh(bid);
    end
    
    %Plot routine
    sp(1) = subplot(211);
    if j == 1
        plot(bdt,zeros(size(bdt)),':k','linewidth',1.5),hold on
    end
    plot(bdt,(bh).*1000,'color',c(j,:),'linewidth',1.5)
    sp(2) = subplot(212);
    p(j) = plot(bdt,eps,'color',c(j,:),'linewidth',1.5); hold on
end
