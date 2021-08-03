%Try plotting the u* ratios by canopy height
clear
load('d:\Mekong_W2015\DataAnalysis\Paper2\RS_estimates\NormFrictionU.mat')
f1 = figure(1);
set(f1,'PaperOrientation','portrait',...
    'position',[400 100 700 600]);
set(gcf,'color','w','paperpositionmode','auto')
symb = {'x';'o';'+';'^'};
c = [0 0 0;1 0 0;0 0 1];
fn = fieldnames(mydata);
hold on
hc = [0.64 0.59 0.61 0.6]; %m, height of canopy
vph = [0.062 0.063 0.061;0.2 0.2 0.2;0.5 0.5 0.5;0.07 0 0];
for i = 1:4
    dfn = fieldnames(mydata.(fn{i}));
    if i < 4
        for ii = 1:3
            vh = vph(i,ii)-0.04-linspace(0.001,0.03,35);
            for j = 1:length(mydata.(fn{i}).(dfn{ii}).time)
                id = mydata.(fn{i}).(dfn{ii}).uid(j);
                zhc = vh(id)/hc(i);
                plot(mydata.(fn{i}).(dfn{ii}).usnorm(j),zhc,...
                    symb{i},'color',c(ii,:))
            end
        end
    else
        vh = vph(i,1)-0.04-linspace(0.001,0.03,35);
        for j = 1:length(mydata.(fn{i}).time)
            id = mydata.(fn{i}).uid(j);
            zhc = vh(id)/hc(i);
            plot(mydata.(fn{i}).usnorm(j),zhc,...
                symb{i},'color',c(ii,:))
        end
    end
end
set(gca,'yscale','log')
