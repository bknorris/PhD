%Try plotting the u,v and w ratios against normalized canopy height
clear
load('d:\Mekong_W2015\DataAnalysis\Paper2\AveragedVelocities\NormalizedVels.mat')
f1 = figure(1);
set(f1,'PaperOrientation','portrait',...
    'position',[400 100   1000   400]);
set(gcf,'color','w','paperpositionmode','auto')
symb = {'o','d','p'};symb = repmat(symb,3,1);
sy = repmat({'^'},1,3);symb = [symb;sy];
c = [0 0 0;1 0 0;0 0 1];
fn = fieldnames(mydata);
hold on
hc = [0.64 0.59 0.61 0.6]; %m, height of canopy
vph = [0.062 0.063 0.061;0.2 0.2 0.2;0.5 0.5 0.5;0.07 0.42 0];
sp(1) = subplot(131);
for i = 1:4
    dfn = fieldnames(mydata.(fn{i}));
    for ii = 1:length(dfn)
        vh = vph(i,ii)-0.04-linspace(0.001,0.03,35);
        for j = 1:length(mydata.(fn{i}).(dfn{ii}).urav)
            id = mydata.(fn{i}).(dfn{ii}).uiav(j);
            zhc = vh(id)/hc(i);
            plot(mydata.(fn{i}).(dfn{ii}).urav(j),zhc,...
                symb{i,ii},'color',c(ii,:)),hold on
        end
    end
end
hold off;
sp(2) = subplot(132);
for i = 1:4
    dfn = fieldnames(mydata.(fn{i}));
    for ii = 1:length(dfn)
        vh = vph(i,ii)-0.04-linspace(0.001,0.03,35);
       for j = 1:length(mydata.(fn{i}).(dfn{ii}).vrav)
            id = mydata.(fn{i}).(dfn{ii}).viav(j);
            zhc = vh(id)/hc(i);
            plot(mydata.(fn{i}).(dfn{ii}).vrav(j),zhc,...
                symb{i,ii},'color',c(ii,:)),hold on
        end
    end
end
hold off
sp(3) = subplot(133);
for i = 1:4
    dfn = fieldnames(mydata.(fn{i}));
    for ii = 1:length(dfn)
        vh = vph(i,ii)-0.04-linspace(0.001,0.03,35);
        for j = 1:length(mydata.(fn{i}).(dfn{ii}).wrav)
            id = mydata.(fn{i}).(dfn{ii}).wiav(j);
            zhc = vh(id)/hc(i);
            plot(mydata.(fn{i}).(dfn{ii}).wrav(j),zhc,...
                symb{i,ii},'color',c(ii,:)),hold on
        end
    end
end
hold off
% set(sp,'yscale','log')
