%plot the pdf of the bottom trace relative to the initial bed elevation
clear,close all
load('e:\Mekong_W2015\DataAnalysis\Paper3\Wavelet\CombEventdata_flood.mat');
data.flood = dat;
load('e:\Mekong_W2015\DataAnalysis\Paper3\Wavelet\CombEventdata_ebb.mat');
data.ebb = dat;clear dat
fn = fieldnames(data);
ffn = fieldnames(data.flood);
f1 = figure(1);
set(f1,'PaperOrientation','portrait',...
    'position',[400 100   800   500],...
    'renderer','painters');
cl = [0.7 0.7 0.7;0.4 0.4 0.4;0.1 0.1 0.1];
symb = {'o','s','^'};
sp = zeros(2,2);
w = [1 2;3 4];
for j = 1:2
    sp(w(j,1)) = subplot(2,2,w(j,1));
    plot(zeros(10,1),linspace(-2,100,10),'-',...
        'linewidth',1,'color','k'),hold on
    for i = 1:3
        bdn = [data.(fn{j}).(ffn{i}).ig.bdmed];
        bins = linspace(-0.1,0.1,20);
        N = hist(bdn,bins);
        p(i) = plot(bins,(N./numel(bdn))*100,...
            '-',...
            'marker',symb{i},...
            'color',cl(i,:),...
            'linewidth',1.5);
    end
    sp(w(j,2)) = subplot(2,2,w(j,2));
    plot(zeros(10,1),linspace(-2,100,10),'-',...
        'linewidth',1,'color','k'),hold on
    for i = 1:3
        deltabd = [data.(fn{j}).(ffn{i}).ig.deltbd];
        bins = linspace(-0.1,0.1,20);
        N = hist(deltabd,bins);
        [h,P,st] = chi2gof(deltabd,'Alpha',0.05);
        if h == 0
            disp([fn{j} ' ' ffn{i} ' does not reject the null at a 5% level'])
        end
        fprintf('Number of events: %0.0f\n',numel(deltabd))
        p(i) = plot(bins,(N./numel(deltabd))*100,...
            '-',...
            'marker',symb{i},...
            'color',cl(i,:),...
            'linewidth',1.5);
    end
end
set(sp(1),'position',[0.11 0.56 0.4 0.35],...
    'ylim',[-2 80],'xlim',[-0.06 0.06],...
    'xticklabel',[])
set(sp(2),'position',[0.56 0.56 0.4 0.35],...
    'ylim',[-2 80],'xlim',[-0.04 0.04],...
    'xticklabel',[],'yticklabel',[])
set(sp(3),'position',[0.11 0.14 0.4 0.35],...
    'ylim',[-2 80],'xlim',[-0.06 0.06])
set(sp(4),'position',[0.56 0.14 0.4 0.35],...
    'ylim',[-2 80],'xlim',[-0.04 0.04],...
    'yticklabel',[])
xlabel(sp(3),'Elevation relative to initial [m]')
xlabel(sp(4),'Net elevation change over event [m]')
ylabel(sp(1),'Percent occurence')
ylabel(sp(3),'Percent occurence')
leg = legend(p,{'Mudflat','Fringe','Forest'});
set(leg,'position',[0.91 0.84 0.005 0.005])
prettyfigures('text',12,'labels',13,'box',1)
sfdir = 'd:\GradSchool\DataAnalysis\Paper3\Figures\';
% export_fig([sfdir 'AllData_bdevel_pdf'],'-pdf')