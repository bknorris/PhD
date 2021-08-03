%Plot figure for Steve, RS versus epsilon and dudz.
clear, close all
load('D:\Projects\Mekong_W2015\DataAnalysis\Paper2\RS_estimates\RStress_10min_adapbdadj.mat')
load('d:\Projects\Mekong_W2015\DataAnalysis\Paper2\EddyViscosity\U_z_gradient_v4.mat')
ddir = 'd:\Projects\Mekong_W2015\DataAnalysis\Paper2\Turbulence\';
files = {'7March2015_TKEbdadj.mat';'8March2015_VelsTKE.mat';'10March2015_VelsTKE.mat';...
    '14March2015a_TKEbdadj.mat'};
dn = {'day1';'day2';'day3';'day4'};
symb = {'o';'o';'^'};
c = [0 0 0;0.3 0.3 0.3;0.6 0.6 0.6;0.9 0.9 0.9];
xs = zeros(4,3);
ys = zeros(4,3);
p = zeros(4,1);
for i = 1:4
    load([ddir files{i}])
    fn = fieldnames(Stat);
    if i == 1 || i == 4
        bin = 1:5;
    else
        bin = 1:30;
    end
    for j = 1:3
        %plot labeling
        if strcmp(dn{i},'day3')
            s = 3;
        elseif strcmp(dn{i},'day4')
            s = j;
        else
            s = 1;
        end
        %Variables to plot
        uw = -1.*nanmean(RS.(dn{i}).(fn{j}).uw(:,bin),2);
        uwt = RS.(dn{i}).(fn{j}).time;
        dudz = slope.(dn{i}).(fn{j}).dudz;
        eps = nanmean((Stat.(fn{j}).z1.E(bin,:)+Stat.(fn{j}).z2.E(bin,:))./2)';
        ept = Stat.(fn{j}).time;
        if i == 2
            uw = uw(50:80);
            uwt = uwt(50:80);
            eps = eps(60:140);
            ept = ept(60:140);
            if j == 1
                dudz = dudz*-1;
            end
        end
        x = (eps.^(1/3)).*dudz;
        newx = zeros(length(uwt),1);
        for k = 1:length(uwt)-1
            id = find(ept>=uwt(k)&ept<=uwt(k+1));
            newx(k) = nanmean(x(id));
        end
        newx(end) = x(end);xstd = nanstd(newx);
        y = uw;ystd = nanstd(y);
        
        p(i) = plot(nanmean(newx)/10,nanmean(y),...
            symb{s},'color','k',...
            'linewidth',1,...
            'markerfacecolor',c(i,:),...
            'markeredgecolor','k',...
            'markersize',6);hold on
        xs(i,j) = nanmean(newx)/10;
        ys(i,j) = nanmean(y);
        
        %         p(i) = plot(newx/10,y,...
        %             symb{s},'color','k',...
        %             'linewidth',1,...
        %             'markerfacecolor',c(i,:),...
        %             'markeredgecolor','k',...
        %             'markersize',6);hold on
        %         eb = ploterr(nanmean(newx),nanmean(y),[ystd ystd],[xstd xstd],...
        %         symb{s},'abshhxy',0.001);hold on
        %         set(eb,...
        %         'color','k',...
        %         'linewidth',1,...
        %         'markerfacecolor',c(i,:),...
        %         'markeredgecolor','k',...
        %         'markersize',6);
    end
end
xs = reshape(xs,12,1);
ys = reshape(ys,12,1);
pf = polyfit(xs,ys,1);
pv = polyval(pf,xs);
p(5) = plot(xs,pv,'-k','linewidth',1.5);
text = sprintf('y = %0.2fx%0.1d',pf(1),pf(2));
leg = legend(p,{'HTA1';'HTA2';'HTA3';'VTA';text},'location','southeast');
xlabel('$[\epsilon^{1/3}]\cdot\delta \overline{u}\slash\delta z$','interpreter','latex')
ylabel('$\overline{-u''w''}$','interpreter','latex')
axis square
prettyfigures('text',12,'labels',14,'box',1)
figdir = 'f:\GradSchool\DataAnalysis\Paper2\WorkingFigures\EddyViscosity\';
export_fig([figdir 'RSvsEddyVisc'],'-pdf','-nocrop')
