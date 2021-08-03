%The multiple linear regression of forcing variables does not identify one
%(or more) variables as the principal variable responsible for the variance
%of the change in bed level. To estimate the relative importance of forcing
%variables, I calculate the principal component regression (PCR) of the
%covariates ('predictors') against the dependent variable (BLE in this
%case).
clear,close all
load('g:\Mekong_W2015\DataAnalysis\Paper3\Wavelet\CombEventData_flood.mat');
% load('g:\Data\Paper3\2015\Wavelet\CombEventData_ebb.mat')
flds = fieldnames(dat);
sdir = 'f:\GradSchool\DataAnalysis\Paper3\WorkingFigures\MLR\';
% cl = brewermap(12,'Paired');
cl = [207 176 126;
    60 166 74;
    4 76 41]./255;
areas = {'Mudflat';'Fringe';'Forest'};
labels = {'Bdmed';'Eventl';'H';'usqd';'taub';'eps'};
explist = {'MA1';'F2F_1';'F2F_2';'F2F2_1';'F2F2_2';'HTA1';'HTA2';...
    'HTA3';'F2F3_1';'F2F3_2';'VTA1';'VTA2'};
symb = {'^';'s';'d';'v';'o';'p';'+';'*';'<';'h';'>';'x'};
pl = zeros(3,1);
figure()
for i = 1:3
    %% PCA - Loading Plot
    exp = dat.(flds{i}).wave.exp;
    eventl = dat.(flds{i}).wave.eventl;
    bdmed = dat.(flds{i}).wave.bdmed;
    Hs = dat.(flds{i}).wave.sigh;
    H = dat.(flds{i}).wave.depth;
    uorb = dat.(flds{i}).wave.orbwv;
    usqd = dat.(flds{i}).wave.usqd;
    taub = dat.(flds{i}).wave.taub;
    eps = dat.(flds{i}).wave.eps;
    data = [bdmed eventl H usqd taub eps];
    %Run PCA on forcing variables
    w = 1./var(data);
    [wcoeff,score,latent,tsquared,explained] = pca(data,...
        'VariableWeights','variance');
    coefforth = inv(diag(std(data)))*wcoeff;
    %Print results of PCA - identify variables with W > 0.5
    varids = [find(abs(coefforth(:,1))>0.5);find(abs(coefforth(:,2))>0.5);];
    sprintf('PCs 1 & 2 explain %0.0f percent of the variance',explained(1)+explained(2))
    for j = 1:length(varids)
        disp(['Significant variables: ' labels{varids(j)}])
    end
    %Plot routine
    plot(zeros(10,1),linspace(-1,1,10),'--k'),hold on
    plot(linspace(-1,1,10),zeros(10,1),'--k')
    for j = 1:length(coefforth)
        pl(i)=plot(coefforth(j,1),coefforth(j,2),symb{j},...
            'markersize',8,...
            'color',cl(i,:),...
            'markerfacecolor',cl(i,:));
    end
    text(coefforth(:,1)+0.02,coefforth(:,2),labels)
end
legend(pl,areas,'location','southwest')
axis([-0.6 0.6 -0.8 0.8])
set(gca,'xtick',-0.6:0.2:0.6,...
    'ytick',-0.6:0.2:0.6)

%% The code below plots the scores & loadings for an individual location
%     figure(i)
%     plot(zeros(10,1),linspace(-1,1,10),'--k'),hold on
%     plot(linspace(-1,1,10),zeros(10,1),'--k')
%     expnum = unique(exp);
%     pl = size(expnum);
%     for j = 1:length(expnum)
%         pl(j) = plot(score(exp==expnum(j),1)./10,score(exp==expnum(j),2)./10,...
%             symb{expnum(j)},...
%             'color',cl(expnum(j),:),...
%             'markersize',6,...
%             'markerfacecolor',cl(expnum(j),:));
%     end
%     plot(coefforth(:,1),coefforth(:,2),'ok',...
%         'markersize',6,...
%         'markerfacecolor','k')
%     text(coefforth(:,1)+0.02,coefforth(:,2),labels)
%     legend(pl,explist{expnum})
%     if i == 2
%         axis([-0.6 0.8 -1 1])
%     else
%         axis([-0.6 0.8 -0.8 0.8])
%     end
%     xlabel('PC1')
%     ylabel('PC2')
%     title(['PCA - ' flds{i}])
%     export_fig([sdir 'PCA_Forest'],'-png')
%% PCR to regress bed level against important components
% exp = dat.(flds{i}).wave.exp;
% bdmed = dat.(flds{i}).wave.bdmed;
% Hs = dat.(flds{i}).wave.sigh;
% uorb = dat.(flds{i}).wave.orbwv;
% eps = dat.(flds{i}).wave.eps;
% data = [Hs uorb eps];
% w = 1./var(data);
% [wcoeff,score,latent,tsquared,explained] = pca(data,...
%     'VariableWeights',w);
% coefforth = inv(diag(std(data)))*wcoeff;
% betaPCR = regress(bdmed-mean(bdmed),score(:,1:2));
% betaPCR = coefforth(:,1:2)*betaPCR;
% betaPCR = [mean(bdmed) - mean(data)*betaPCR; betaPCR];
% yfitPCR = [ones(length(data),1) data]*betaPCR;
% TSS = sum((bdmed-mean(bdmed)).^2);
% RSS_PCR = sum((bdmed-yfitPCR).^2);
% rsquaredPCR = 1 - RSS_PCR/TSS;
% figure(i)
% plot(bdmed,yfitPCR,'^r')
%% Simple Auto-correlation
%      figure(i)
%      C = corr(data,data);
%     imagesc(C)
%     set(gca,'xticklabel',labels,'yticklabels',labels)
%     caxis([-1 1]),cb = colorbar;
%     ylabel(cb,'Correlation')
%     title(['Auto Correlation of ' flds{i}])
%     colormap(cl)
%     prettyfigures('text',12,'labels',13,'box',1)
%     export_fig([sdir 'CorrCoef_' flds{i}],'-png')

