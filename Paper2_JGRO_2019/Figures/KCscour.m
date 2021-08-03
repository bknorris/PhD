%Plot KC # versus scour depth for the final fig of Paper 2.
clear
load('D:\Projects\Mekong_W2015\DataAnalysis\Paper2\WVStats_HTA_VTA.mat')
wv = data;
load('D:\Projects\Mekong_W2015\DataAnalysis\Paper2\BDtrace_HTA1_2.mat')
bd = data;
load('D:\Projects\Mekong_W2015\DataAnalysis\Paper2\HTAVTA_UwUcStLt.mat')

f1 = figure(1);
set(f1,'PaperOrientation','portrait',...
    'position',[400 100   800   600]);
names = {'vp1';'vp2';'vp3'};
fn = fieldnames(data.day1);
symb = {'o';'d';'p'};
lines = {'-';'--';'-.'};
tstep = datenum(0,0,0,0,30,0);
sp = zeros(4,1);
cl = flipud([0.1 0.1 0.1;0.5 0.5 0.5;0.7 0.7 0.7]);
pl = zeros(3,1);
for i = 1:3
    sp(1) = subplot(121);
    time = bd.day1.(names{i}).time;
    time2 = data.day1.(fn{i}).time;
    lvl = bd.day1.(names{i}).lvl.*-1;
    
    Tr = wv.V1.Tr;Tr(isnan(Tr)) = 0;
    T = interp1(1:length(wv.V1.time2),Tr,1:length(time));
    Uw = interp1(1:length(time2),data.day1.(fn{i}).Uw,1:length(time));
    d = data.day1.(fn{i}).d(1);
    KC = (Uw.*T)./d;
    KC(KC<6) = NaN;
    n = round(length(KC)/20);
    avt = [1 n:n:length(KC) length(KC)];
    for ii = 1:length(avt)-1
        kc = mean(KC(avt(ii):avt(ii+1)));
        Sd = mean(lvl(avt(ii):avt(ii+1)))/d;
        plot(kc,Sd,symb{i},...
            'color',cl(i,:),'linewidth',1.5), hold on
    end
    
    sp(2) = subplot(122);
    time = bd.day2.(names{i}).time;
    time2 = data.day2.(fn{i}).time;
    lvl = bd.day2.(names{i}).lvl.*-1;
    Tr = wv.V2.Tr;Tr(isnan(Tr)) = 0;
    T = interp1(1:length(wv.V2.time2),Tr,1:length(time));
    Uw = interp1(1:length(time2),data.day2.(fn{i}).Uw,1:length(time));
    d = data.day2.(fn{i}).d(1);
    KC = (Uw.*T)./d;
    KC(KC<6) = NaN;
    n = round(length(KC)/20);
    avt = [1 n:n:length(KC) length(KC)];
    for ii = 1:length(avt)-1
        kc = mean(KC(avt(ii):avt(ii+1)));
        Sd = mean(lvl(avt(ii):avt(ii+1)))/d;
        plot(kc,Sd,symb{i},...
            'color',cl(i,:),'linewidth',1.5), hold on
    end
end

