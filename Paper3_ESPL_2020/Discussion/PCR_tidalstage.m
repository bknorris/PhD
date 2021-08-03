clear,close all
load('d:\Mekong_W2015\DataAnalysis\Paper3\CmbData\CombAllData_flood.mat');
data.flood = dat;
load('d:\Mekong_W2015\DataAnalysis\Paper3\CmbData\CombAllData_ebb.mat');
data.ebb = dat;clear dat
fn = fieldnames(data);
sfdir = 'e:\GradSchool\DataAnalysis\Paper3\Figures\';
flds = {'mud';'fringe';'forest'};
stage = {'Low';'Medium-Low';'Medium-High';'High'};
labels = {'Hs';'u-squared';'ubr';'taub';'eps'};
for i = 1:2
    for ii = 1:3
        %% PCA - compute PC1 & PC2
        time = data.(fn{i}).(flds{ii}).time;
        bdist = data.(fn{i}).(flds{ii}).bdist;
        Hs = data.(fn{i}).(flds{ii}).Hs;
        H = data.(fn{i}).(flds{ii}).depth;
        ubr = data.(fn{i}).(flds{ii}).ubr;
        usqd = data.(fn{i}).(flds{ii}).usqd;
        taub = data.(fn{i}).(flds{ii}).tmax;
        eps = data.(fn{i}).(flds{ii}).eps;
        %Parse model data by water depth
        dbins = [0 0.4 0.6 0.8 1 1.4];
        for k = 1:4
            fprintf([stage{k} ' tidal stage\n***\n'])
            idx = find(H>=dbins(k)&H<=dbins(k+1));
            model = [Hs(idx); usqd(idx); ubr(idx); taub(idx); eps(idx)]';
            output = pca_model(model,2,'auto');
            varids = [find(abs(output.L(:,1))>0.5);find(abs(output.L(:,2))>0.5);];
            sprintf('PCs 1 & 2 explain %0.0f percent of the variance',(output.exp_var(1)+output.exp_var(2))*100)
            for j = 1:length(varids)
                disp(['Significant variables: ' labels{varids(j)}])
            end
            figure(k)
            bar(output.L)
            set(gca,'xticklabel',labels)
            title([stage{k} ' ' upper(fn{i}) ' Tide - ' upper(flds{i})])
            if k == 4
                tilefigs()
                disp('Press any key to continue')
                pause
                close all
            end
        end
    end
end

