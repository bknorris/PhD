% Plot depth profiles of removed RS estimates for the three instruments of
% HTA1-3 and VTA.
clear
% load('D:\Mekong_W2015\DataAnalysis\Paper2\RStress_10min_adapbdadj.mat')
% pct = [33 46 69;51 35 45;41 52 45;36 55 0];
% fn = fieldnames(RS);
% for k = 1:4
%     dfn = fieldnames(RS.(fn{k}));
%     nanprof = zeros(35,1);
%     rng('default')
%     for i = 1:3
%         uw = RS.(fn{k}).(dfn{i}).uw;
%         npts = numel(uw);
%         nnan = (pct(k,i)/100)*npts;
%         cnan = length(find(isnan(uw)));
%         rnan = nnan-cnan;
%         [nsamp,~] = size(uw);
%         nancol = round(rnan/nsamp);
%         pctnan = length(find(isnan(uw)))/numel(uw);
%         ct = 1;
%         while pctnan < pct(k,i)/100
%             row = randi(nsamp,1,1);
%             if mod(ct,2)
%                 rdbeg = round(sqrt(35)*randn(1,35) + 1);
%                 rdend = round(sqrt(35)*randn(1,35) + 35);
%                 ind = [rdbeg(rdbeg>0) rdend(rdend<=35)]; 
%             else
%                 ind = randi(35,1,25);
%             end
%             uw(row,ind) = NaN;
%             pctnan = length(find(isnan(uw)))/numel(uw);
%             ct = ct+1;
%         end
%         RS.(fn{k}).(dfn{i}).uw = uw;
%         for ii = 1:35
%             nanprof(ii) = length(find(isnan(uw(:,ii))))/numel(uw(:,ii));
%             percent nans
%         end
%         RS.(fn{k}).(dfn{i}).pctnan = nanprof;
%         
%     end
% end
% sfd = 'd:\Mekong_W2015\DataAnalysis\Paper2\RS_estimates\Final\';
% save([sfd 'RStress_10min_bdadj_f.mat'],'RS','-v7.3')
load('d:\Mekong_W2015\DataAnalysis\Paper2\RS_estimates\Final\RStress_10min_bdadj_f.mat')
savefigdir = 'c:\Users\Bnorr\Documents\GradSchool\Writing\JGR_2018_Norris\Figures\Draft2_Figures\Versions\V2\';

%Plot Routine
f1 = figure(1);
set(f1,'PaperOrientation','portrait',...
    'position',[400 100   1000   400]);
set(gcf,'color','w','paperpositionmode','auto')
sp = zeros(4,1);
dn = fieldnames(RS);
% hab = [0.062 0.063 0.061;
%     0.242 0.240 0.240;
%     0.550 0.550 0.550;
%     0.07 0.416 0.806];
% hc = [0.58 0.58 0.58 0.6];
symb = {'o';'d';'p'};
line = {'-';'--';'-.'};
titles = {'HTA1';'HTA2';'HTA3';'VTA'};
cl = flipud([0.1 0.1 0.1;0.5 0.5 0.5;0.7 0.7 0.7]);
for i = 1:4
    fn = fieldnames(RS.(dn{i}));
    sp(i) = subplot(1,4,i);
    for j = 1:3
        pctnan = RS.(dn{i}).(fn{j}).pctnan;
        %         z = hab(i,j)-0.04-linspace(0,0.03,35);
        %             zhc = z./hc(i);
        z = (0.001:0.001:0.035)+0.04;
        markx = (1-pctnan(1:5:end))*100;
        marky = z(1:5:end);
        plot((1-pctnan)*100,z,line{j},...
            'color',cl(j,:),...
            'linewidth',1.5), hold on
        plot(markx,marky,symb{j},...
            'linewidth',1.5,...
            'markerfacecolor',cl(j,:),...
            'markeredgecolor','k',...
            'markersize',6)
        grid on
    end
    xlabel('% of bursts')
    title(titles{i})
end
hold off
%global plot adjustments
set(sp,'xlim',[0 100],...
    'xtick',0:25:100,...
    'ydir','reverse',...
    'ylim',[0.04 0.075],...
    'gridlinestyle',':')
set([sp(2) sp(3) sp(4)],...
    'yticklabel',[])
ylabel(sp(1),'Distance from sensor (m)')
%plot positioning
set(sp(1),'position',[0.1 0.18 0.19 0.72])
set(sp(2),'position',[0.33 0.18 0.19 0.72])
set(sp(3),'position',[0.56 0.18 0.19 0.72])
set(sp(4),'position',[0.79 0.18 0.19 0.72])
prettyfigures('text',13,'labels',14,'box',1,'tlength',[0.025 0.025])
export_fig([savefigdir 'RSremovedPts'],'-pdf')

