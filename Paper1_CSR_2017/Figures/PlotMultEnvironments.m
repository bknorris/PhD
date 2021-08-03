%Try dividing up the dataset into mudflat, fringe and forest, and plotting
%variables to look for patterns.
clear
BigDataConglomerator
plotsmscale = 1;
plotlgscale = 0;

%define boundaries
mudflat = [-80 -20]; fringe = [-20 20]; forest = [20 100];
savefigdir = 'd:\Projects\Mekong_W2015\Figures\Paper1\';

if plotsmscale
    id1 = find(X >= mudflat(1) & X <= mudflat(2));
    id2 = find(X >= fringe(1) & X <= fringe(2));
    id3 = find(X >= forest(1) & X <= forest(2));
    
    %boxplot phi values on/off the mudflat
    mud = eps(id1);fr = eps(id2);fo = eps(id3);
    bpx = NaN(length(fr),3);
    bpx(1:length(mud),1) = mud;
    bpx(1:length(fr),2) = fr;
    bpx(1:length(fo),3) = fo;
    bpx(bpx == 0) = NaN;
    
    %colors (mudflat,fringe,forest)
    c = [137 67 54;107 206 82;65 134 78]./255;
    %Figure Plotting Routine
    f1 = figure(1);
    set(f1,'PaperOrientation','portrait',...
        'position',[400 200   600   500]);
    set(gcf,'color','w','PaperPositionMode','auto')
    bp = boxplot(bpx,'colors','k','Labels',{'Mudflat','Fringe','Forest'});
    set(bp,'linewidth',2);
    set(gca,'Yscale','log','Ylim',[1E-7 1E-1],...
        'LineWidth',1.5,'FontSize',12,'FontName','Cambria')
    ylabel('\epsilon (WKg^-^1)','FontSize',12,'FontName','Cambria')
    set(findobj(gca,'Type','text'),'FontSize',12,'FontName','Cambria')
    % export_fig([savefigdir 'MultEnviBoxPlot'],'-jpeg','-nocrop')
    
    %try plotting phi/epsilon for the three locations
    f2 = figure(2);
    set(f2,'PaperOrientation','portrait',...
        'position',[400 200   1000   500]);
    set(gcf,'color','w','PaperPositionMode','auto')
    p(1) = subplot(131);
    plot(phi(id1),mud,'d','Color',c(1,:),...
        'markersize',12,'MarkerFaceColor','w','LineWidth',1.5)
    xlabel('\phi','FontSize',14,'FontName','Cambria')
    ylabel('\epsilon (WKg^-^1)','FontSize',14,'FontName','Cambria')
    title('Mudflat','FontSize',14,'FontName','Cambria')
    p(2) = subplot(132);
    %linreg of phi and eps
    frphi = phi(id2);frphi = frphi(~isnan(fr));fr = fr(~isnan(fr));
    pf = polyfit(frphi,log(fr),1);
    x1= linspace(-0.05,0.1,length(frphi));
    y1 = exp(pf(1)*x1+pf(2));
    plot(x1,y1,'-k','LineWidth',1.5), hold on
    plot(frphi,fr,'d','Color',c(2,:),...
        'markersize',12,'MarkerFaceColor','w','LineWidth',1.5)
    xlabel('\phi','FontSize',14,'FontName','Cambria')
    title('Fringe','FontSize',14,'FontName','Cambria')
    p(3) = subplot(133);
    %linreg of phi and eps
    fophi = phi(id3);fophi = fophi(~isnan(fo));
    fo = fo(~isnan(fo));fo = fo(fo~=0);fophi = fophi(fo~=0);
    pf = polyfit(fophi,log(fo),1);
    x1= linspace(-0.05,0.1,length(fophi));
    y1 = exp(pf(1)*x1+pf(2));
    plot(x1,y1,'-k','LineWidth',1.5), hold on
    plot(fophi,fo,'d','Color',c(3,:),...
        'markersize',12,'MarkerFaceColor','w','LineWidth',1.5)
    xlabel('\phi','FontSize',14,'FontName','Cambria')
    title('Forest','FontSize',14,'FontName','Cambria')
    set(p,'XLim',[-0.01 0.1],'XTick',0:0.05:0.1,...
        'YLim',[1E-7 1E-1],'Yscale','log',...
        'LineWidth',1.5,'FontSize',14,'FontName','Cambria')
    set([p(2) p(3)],'YTickLabel',[])
    %positioning
    set(p(1),'position',[0.1 0.12 0.26 0.8])
    set(p(2),'position',[0.39 0.12 0.26 0.8])
    set(p(3),'position',[0.68 0.12 0.26 0.8])
    % export_fig([savefigdir 'MultEnviPhiEps'],'-jpeg','-nocrop')
    
    %try plotting significant wave height/epsilon for the three locations
    f3 = figure(3);
    set(f3,'PaperOrientation','portrait',...
        'position',[400 200   1000   500]);
    set(gcf,'color','w','PaperPositionMode','auto')
    p(1) = subplot(131);
    mudhs = Hs(id1);mudhs = mudhs(~isnan(mud));mud = mud(~isnan(mud));
    mud(isnan(mudhs)) = [];mudhs(isnan(mudhs)) = [];
    pf = polyfit(mudhs,log(mud),1);
    x1= linspace(0,0.8,length(mudhs));
    y1 = exp(pf(1)*x1+pf(2));
    plot(x1,y1,'-k','LineWidth',1.5), hold on
    plot(mudhs,mud,'o','Color',c(1,:),...
        'markersize',8,'MarkerFaceColor','w','LineWidth',1.5)
    xlabel('H_s (m)','FontSize',14,'FontName','Cambria')
    ylabel('\epsilon (WKg^-^1)','FontSize',14,'FontName','Cambria')
    title('Mudflat','FontSize',14,'FontName','Cambria')
    p(2) = subplot(132);
    %linreg of phi and eps
    frhs = Hs(id2);frhs = frhs(~isnan(fr));fr = fr(~isnan(fr));
    fr(isnan(frhs)) = [];frhs(isnan(frhs)) = [];
    pf = polyfit(frhs,log(fr),1);
    x1= linspace(0,0.8,length(frhs));
    y1 = exp(pf(1)*x1+pf(2));
    plot(x1,y1,'-k','LineWidth',1.5), hold on
    plot(frhs,fr,'o','Color',c(2,:),...
        'markersize',8,'MarkerFaceColor','w','LineWidth',1.5)
    xlabel('H_s (m)','FontSize',14,'FontName','Cambria')
    title('Fringe','FontSize',14,'FontName','Cambria')
    p(3) = subplot(133);
    %linreg of phi and eps
    fohs = Hs(id3);fohs = fohs(~isnan(fo));
    fo = fo(~isnan(fo));fo = fo(fo~=0);fohs = fohs(fo~=0);
    pf = polyfit(fohs,log(fo),1);
    x1= linspace(0,0.8,length(fophi));
    y1 = exp(pf(1)*x1+pf(2));
    plot(x1,y1,'-k','LineWidth',1.5), hold on
    plot(fophi,fo,'o','Color',c(3,:),...
        'markersize',8,'MarkerFaceColor','w','LineWidth',1.5)
    xlabel('H_s (m)','FontSize',14,'FontName','Cambria')
    title('Forest','FontSize',14,'FontName','Cambria')
    set(p,'XLim',[0 0.8],'XTick',0:0.2:0.8,...
        'YLim',[1E-7 1E-1],'Yscale','log',...
        'LineWidth',1.5,'FontSize',14,'FontName','Cambria')
    set([p(2) p(3)],'YTickLabel',[])
    %positioning
    set(p(1),'position',[0.1 0.13 0.26 0.8])
    set(p(2),'position',[0.39 0.13 0.26 0.8])
    set(p(3),'position',[0.68 0.13 0.26 0.8])
    % export_fig([savefigdir 'MultEnviHsEps'],'-jpeg','-nocrop')
end

%try averaging phi over larger distances of say 10m (and average eps to
%make the arrays the same size).
if plotlgscale
    Xs = [reshape(X,153,1);xqq];Phi = [reshape(phi,153,1); phiqq];
    muddist = abs(mudflat(2)-mudflat(1))/3;
    frdist = abs(fringe(2)-fringe(1))/4;
    fodist = abs(forest(2)-forest(1))/2;
    %select and average the data
    mud = mudflat(1):muddist:mudflat(2);
    mphiavg = zeros(3,1); %phi
    mepsavg = zeros(3,1); %epsilon
    mstd = zeros(3,1); %stdev
    for i = 1:length(mud)-1
        id = find(Xs >= mud(i) & Xs <= mud(i+1));
        mphiavg(i) = nanmean(Phi(id)); %#ok<*FNDSB>
        id = find(X >= mud(i) & X <= mud(i+1));
        mepsavg(i) = nanmean(eps(id));
        mstd(i) = nanstd(eps(id));
    end
    fr = fringe(1):frdist:fringe(2);
    frphiavg = zeros(4,1); %phi
    frepsavg = zeros(4,1); %epsilon
    frstd = zeros(3,1); %stdev
    for i = 1:length(fr)-1
        id = find(Xs >= fr(i) & Xs <= fr(i+1));
        frphiavg(i) = nanmean(Phi(id));
        id = find(X >= fr(i) & X <= fr(i+1));
        frepsavg(i) = nanmean(eps(id));
        frstd(i) = nanstd(eps(id));
    end
    fo = forest(1):fodist:forest(2);
    fophiavg = zeros(2,1); %phi
    foepsavg = zeros(2,1); %epsilon
    fostd = zeros(2,1); %stdev
    for i = 1:length(fo)-1
        id = find(Xs >= fo(i) & Xs <= fo(i+1));
        if i == 2
            fophiavg(i) = 0.015;
        else
            fophiavg(i) = nanmean(Phi(id));
        end
        id = find(X >= fo(i) & X <= fo(i+1));
        foepsavg(i) = nanmean(eps(id));
        fostd(i) = nanstd(eps(id));
    end
    
    f1 = figure(1);
    set(f1,'PaperOrientation','portrait',...
        'position',[400 200   1000   500]);
    set(gcf,'color','w','PaperPositionMode','auto')
    Xd = [-70 -50 -30 -15 -5 5 15 40 80]';
    totalphi = [mphiavg;frphiavg;fophiavg];
    totaleps = [mepsavg;frepsavg;foepsavg];
    totalstd = [mstd;frstd;fostd];
    p(1) = plot(Xd,totalphi,'o','MarkerEdgeColor','k',...
        'MarkerFaceColor','w','MarkerSize',10,...
        'LineWidth',1.5);hold on
    %vertical line at zero
    yl = linspace(-1,1,100);xl = zeros(1,length(yl));
    p(2) = plot(xl,yl,'--k','LineWidth',1.5);
    axis([-80 100 0 0.03])
    ax1 = gca;
    ax1_pos = get(ax1,'position'); % position of first axes
    ax2 = axes('Position',ax1_pos,...
        'XAxisLocation','top',...
        'YAxisLocation','right',...
        'Color','none');
    hold on
    errorbar(Xd,totaleps,totalstd,'.',...
        'Color','r','MarkerSize',8,'LineWidth',1.5,...
        'parent',ax2)
    set(ax2,'YLim',[-1.5E-3 5E-3],'xticklabel',[],'LineWidth',1.5,'box','on',...
        'FontSize',14,'FontName','Cambria')
    ylabel(ax2,'\epsilon (WKg^-^1)','FontSize',14,'FontName','Cambria')
    xlabel(ax1,'Cross-Shore Distance (m)','FontSize',14,'FontName','Cambria')
    ylabel(ax1,'\phi','FontSize',14,'FontName','Cambria')
    set(ax1,'FontSize',14,'FontName','Cambria')
    %    errorbar(Xd,totaleps,totalstd,'.',...
    %         'Color','k','MarkerSize',8,'LineWidth',1.5)
    %     %vertical line at zero
    yl = linspace(-1,1,100);xl = zeros(1,length(yl));
    p(2) = plot(xl,yl,'--k','LineWidth',1.5);
    axis([-80 100 -1.5E-3 4E-3])
    set(gca,'YTick',-1.5E-3:1E-3:4E-3,'XTick',-80:20:100,...
        'LineWidth',1.5,'FontSize',14,'FontName','Cambria')
    export_fig([savefigdir 'AvgPhiEps'],'-jpeg','-nocrop')

end