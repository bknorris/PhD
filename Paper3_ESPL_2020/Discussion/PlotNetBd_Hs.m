clear,close all
load('f:\Mekong_W2015\DataAnalysis\Paper3\Wavelet\CombEventdata_flood.mat');
data.flood = dat;
load('f:\Mekong_W2015\DataAnalysis\Paper3\Wavelet\CombEventdata_ebb.mat');
data.ebb = dat;clear dat
fn = fieldnames(data);
sfdir = 'g:\GradSchool\DataAnalysis\Paper3\Figures\';
flds = {'mud';'fringe';'forest'};
band = {'wave';'ig'};
cl = [207 176 126;
    60 166 74;
    4 76 41]./255;
symb = {'o';'d';};
figure(1),hold on
c = 1;
for i = 1:2
    disp(fn{i})
    for j = 1:length(flds)
        subplot(2,3,c)
        deltbd = abs([data.(fn{i}).(flds{j}).wave.deltbd; data.(fn{i}).(flds{j}).ig.deltbd]);
        deltbd(deltbd>0.1)=0; %remove pesky lg. datapoints
        Hs = [data.(fn{i}).(flds{j}).wave.sigh; data.(fn{i}).(flds{j}).ig.sigh];
        deltbd(Hs<0.1) = [];
        Hs(Hs<0.1) = [];
        CI = (std(deltbd')./sqrt(size(deltbd,2))).*1.96;
        id = find(deltbd>=CI); %these should be the outliers
        f = fit(deltbd(id),Hs(id),'poly2');
        plot(f,deltbd(id),Hs(id))
        c = c+1;
        set(gca,'ylim',[0 0.6],'xlim',[0 0.1])
    end
end
suplabel('Net Elevation Change over Event (m)','x')
suplabel('H_s (m)','y')

figure(2),hold on
c = 1;
for i = 1:2
    disp(fn{i})
    for j = 1:length(flds)
        subplot(2,3,c)
        deltbd = abs([data.(fn{i}).(flds{j}).wave.deltbd; data.(fn{i}).(flds{j}).ig.deltbd]);
        deltbd(deltbd>0.1)=0; %remove pesky lg. datapoints
        taub = [data.(fn{i}).(flds{j}).wave.taub; data.(fn{i}).(flds{j}).ig.taub];
        deltbd(taub<0.1) = [];
        taub(taub<0.1) = [];
        CI = (std(deltbd')./sqrt(size(deltbd,2))).*1.96;
        id = find(deltbd>=CI); %these should be the outliers
        f = fit(deltbd(id),taub(id),'poly1');
        plot(f,deltbd(id),taub(id))
        c = c+1;
        set(gca,'ylim',[0 1.6],'xlim',[0 0.1])
    end
end
suplabel('Net Elevation Change over Event (m)','x')
suplabel('H_s (m)','y')
prettyfigures