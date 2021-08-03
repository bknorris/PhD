%plot the pdf of the bottom trace relative to the initial bed elevation
clear,close all
load('e:\Mekong_W2015\DataAnalysis\Paper3\Wavelet\CombEventdata_flood_v2.mat');
data.flood = dat;
load('e:\Mekong_W2015\DataAnalysis\Paper3\Wavelet\CombEventdata_ebb_v2.mat');
data.ebb = dat;clear dat
fn = fieldnames(data);
ffn = fieldnames(data.flood);
f1 = figure(1);
set(f1,'PaperOrientation','portrait',...
    'position',[400 100   500   600],...
    'renderer','painters');
cl = [0.7 0.7 0.7;0.4 0.4 0.4;0.1 0.1 0.1];
symb = {'o','s','^'};
sp = zeros(2,2);
for j = 1:2
    sp(j) = subplot(2,1,j);
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
set(sp(1),'position',[0.12 0.56 0.8 0.38],...
    'ylim',[-2 80],'xlim',[-0.04 0.04],...
    'xticklabel',[])
set(sp(2),'position',[0.12 0.12 0.8 0.38],...
    'ylim',[-2 80],'xlim',[-0.04 0.04],...
    'xticklabel',{'-40','-20','0','20','40'})
xlabel(sp(2),'Net elevation change over event [mm]')
ylabel(sp(1),'Percent occurence')
ylabel(sp(2),'Percent occurence')
leg = legend(p,{'Mudflat','Fringe','Forest'});
set(leg,'position',[0.84 0.87 0.02 0.02])
prettyfigures('text',12,'labels',13,'box',1)
sfdir = 'd:\GradSchool\DataAnalysis\Paper3\Figures\';
% export_fig([sfdir 'AllData_bdevel_pdf_v2'],'-pdf')