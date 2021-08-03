%estimate total time between wave and IG events
clear,close all
load('f:\Mekong_W2015\DataAnalysis\Paper3\Wavelet\CombEventdata_flood.mat');
data.flood = dat;
load('f:\Mekong_W2015\DataAnalysis\Paper3\Wavelet\CombEventdata_ebb.mat');
data.ebb = dat;clear dat
tide = {'flood';'ebb'};
times = [204+140+105; 80+140+131];
norm = [2 3 2]
for i = 1:2
    area = fieldnames(data.(tide{i}));
    for ii = 1:3
        wevents = data.(tide{i}).(area{ii}).wave.eventl;
        ievents = data.(tide{i}).(area{ii}).ig.eventl;
        disp(['Field area: ' area{ii}])
        fprintf('Total time of wave events: %0.2f minutes\n',(sum(wevents)/60))
        fprintf('Total time of ig events: %0.2f hours\n\n',(sum(ievents)/60/60))
        
        fprintf('Percent of wave events over all experiments: %0.2f %%\n',...
            ((sum(wevents)/60)/times(i))/norm(ii)*100)
        fprintf('Percent of ig events over all experiments: %0.2f %%\n\n',...
            ((sum(ievents)/60)/times(i))/norm(ii)*100)
    end
end