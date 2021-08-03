%MLR using zscore to normalize input values
clear, close all
load('d:\Mekong_W2015\DataAnalysis\Paper3\Wavelet\CombEventData_flood.mat');
flds = fieldnames(dat);
vals = {'H';'Hs';'uorb';'ucub';'taub';'eps';'bdmed'};
stage = {'LL';'ML';'MH';'HH'};c = jet(3);
figno = [3 4 5;6 7 8];
for i = 1 %:length(flds)
    disp(['Running MLR on ' flds{i}])
    exp = [dat.(flds{i}).wave.exp; dat.(flds{i}).ig.exp];
    eventl = [dat.(flds{i}).wave.eventl; dat.(flds{i}).ig.eventl];
    bdmed = [dat.(flds{i}).wave.bdmed; dat.(flds{i}).ig.bdmed];
    Hs = [dat.(flds{i}).wave.sigh; dat.(flds{i}).ig.sigh];
    H = [dat.(flds{i}).wave.depth; dat.(flds{i}).ig.depth];
    uorb = [dat.(flds{i}).wave.orbwv; dat.(flds{i}).ig.orbwv];
    ucub = [dat.(flds{i}).wave.ucub; dat.(flds{i}).ig.ucub];
    taub = [dat.(flds{i}).wave.taub; dat.(flds{i}).ig.taub];
    eps = [dat.(flds{i}).wave.eps; dat.(flds{i}).ig.eps];
%     data = [H Hs uorb ucub taub eps];
    vars = vals([1 4 5 6 7]);
    data = [ones(size(H)) H ucub taub eps];
    y = bdmed;
    [b,bint,r,rint,stats]  = regress(zscore(y),[data(:,1) zscore(data(:,2:end))],0.01);
    fprintf('Variables included: %s\n',vars{:})
    fprintf('R-squared: %0.3f\n',stats(1))
    lma = stepwiselm(zscore(data(:,2:end)),zscore(y),'constant','Upper','linear',...
                    'penter',0.05,'premove',0.1,'varnames',vars)
%     [~,numvar] = size(data);
%     numvar = 2;
%     figure(1)
%     [XL,yl,XS,YS,beta,PCTVAR,mse,stats] = plsregress(zscore(data),zscore(y),numvar);
%     plot(1:numvar,cumsum(100*PCTVAR(2,:)),'-o','color',c(i,:));hold on
%     yfit = [ones(size(data,1),1) data]*beta;
%     figure(2)
%     plot(y,yfit,'o','color',c(i,:)),hold on
% %     TSS = sum((y-mean(y)).^2);
% %     RSS = sum((y-yfit).^2);
% %     Rsquared = 1 - RSS/TSS;
% %     fprintf('R-squared: %0.2f\n',Rsquared)
%     figure(figno(1,i))
%     plot(1:numvar,stats.W,'o-');
%     legend({'c1';'c2';'c3';'c4';'c5'},'Location','NW')
%     set(gca,'xtick',1:1:5,'xlim',[0 5],'xticklabel',vars)
%     xlabel('Predictor');
%     ylabel('Weight');
%     figure(figno(2,i))
%     [axes,h1,h2] = plotyy(0:numvar,mse(1,:),0:numvar,mse(2,:));
%     set(h1,'Marker','o')
%     set(h2,'Marker','o')
%     legend('MSE Predictors','MSE Response')
%     xlabel('Number of Components')
    %     if i < 3
    %         disp('Press any key to continue')
    %         pause
    %     end
    %         data = [bdmed eventl H Hs usqd taub eps];
    %         lma = stepwiselm(data(:,2:end),data(:,1),'constant','Upper','linear',...
    %             'penter',0.05,'premove',0.1,'varnames',vars)
    %         if i < 3
    %             disp('Press any key to continue')
    %             pause
    %         end
end

