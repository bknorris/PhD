%Calculate Depth profiles of mean velocities/Reyolds Stress for the
%Vectrino profilers, from VTA

clear
datdir = 'c:\Users\bkn5\Projects\Mekong_W2015\DataAnalysis\Paper1\';
fname = 'VTAvelocities.mat';
name = 'VTA';
load([datdir fname])
load([datdir 'VTARS.mat'])
fn = fieldnames(VTA);
ph = 42.51; %height of canopy in cm
heading = 102;
start = datenum(2015,03,14,07,00,00);stop = datenum(2015,03,14,11,00,00);tstep = datenum(0,0,0,0,30,0);
hab = [0.07 0.416 0.806]*100;
vname = {'VP1';'VP2';'VP3'};

count = 1;
for i = 2:4
    %velocities
    idx = find(VTA.(fn{i}).time >= start & VTA.(fn{i}).time <= stop);
    u = VTA.(fn{i}).y(idx,:); %along-shore
    v = VTA.(fn{i}).x(idx,:); %cross-shore
    w = (VTA.(fn{i}).z1(idx,:)+VTA.(fn{i}).z2(idx,:))./2;
    t = VTA.(fn{i}).time(idx,:);
    
    %rotate to the cross-shore direction
    rot = heading*pi/180;
    u = u.*(ones(size(u))*cos(rot)) + ...
        v.*(ones(size(v))*sin(rot));
    v = -u.*(ones(size(u))*sin(rot)) + ...
        v.*(ones(size(v))*cos(rot));
    
    h = hab(count); %inst hab (cm)
    rb = VTA.(fn{i}).rb;dat.(fn{i}).rbcm = h-rb(1:end)*100;
    sr = 50;
    intv = 10; %minutes
    avt = intv*sr*60;
    ind = [1 avt:avt:length(u)];
    
    for ii = 1:length(ind)-1
        Uz = u(ind(ii):ind(ii+1),:);
        dat.(fn{i}).U(ii,:) = nanmean(Uz);
        dat.(fn{i}).T(ii,:) = t(ind(ii));
    end
    count = count+1;
end

%%%%%%%% Plotting Routine
fn = fieldnames(dat);
f1 = figure;
set(f1,'PaperOrientation','portrait',...
    'position',[400 100   800   1000]);
set(gcf,'color','w','PaperPositionMode','auto')

for i = 1:3
    sbp1 = [5,3,1]; %subplot order
    pq(i) = subplot(3,2,sbp1(i));
    mVel = mean(dat.(fn{i}).U);Velstd = std(dat.(fn{i}).U);
    ys = (dat.(fn{i}).rbcm)./ph;
    
    zby = linspace(0,2,length(ys));
    zbx = zeros(length(zby),1);
    plot(zbx,zby,'--k','LineWidth',1.5), hold on
%     h1 = herrorbar(mVel(1:3:end)*100,ys(1:3:end),2*100*Velstd(1:3:end));
%     set(h1,'LineWidth',1.5,'Color','k')
    plot(mVel(1:3:end)*100,ys(1:3:end),'Color','k','LineWidth',2,...
        'marker','^','markersize',8,'markerfacecolor',[1 1 1])
    xl = xlabel('U (cm/s)');
    set(xl,'Interpreter','latex','fontsize',14,'FontName','Cambria')
    yl = ylabel('$$z/h_c$$');
    set(yl,'Interpreter','latex','fontsize',14,'FontName','Cambria')
    grid on
    set(gca,'GridLineStyle',':')

    sbp2 = [6,4,2];
    qp(i) = subplot(3,2,sbp2(i));
    mRS = mean(rss.(fn{i}).uw);RSstd = std(rss.(fn{i}).uw);
    
    plot(zbx,zby,'--k','LineWidth',1.5), hold on
%     h2 = herrorbar(mRS(1:3:end),ys(1:3:end),2*RSstd(1:3:end));
%     set(h2,'LineWidth',1.5,'Color','k')
    plot(mRS(1:3:end),ys(1:3:end),'Color','k','LineWidth',2,...
        'marker','^','markersize',8,'markerfacecolor',[1 1 1])
    xl = xlabel(['$$\rho','\langle','\overline{u''v''}','    \rangle$$ (Pa)']);
    set(xl,'Interpreter','latex','fontsize',14,'FontName','Cambria')
    yl = ylabel('$$z/h_c$$');
    set(yl,'Interpreter','latex','fontsize',14,'FontName','Cambria')
    grid on
    set(gca,'GridLineStyle',':')

end

%legend
% leg = legend(pp,'x = -10cm','x = 10cm','x = 20cm');
% set(leg,'position',[0.82 0.58 0.05 0.05],'LineWidth',1.5)
%plot adjustments
set([pq qp],'LineWidth',1.5,'FontSize',14,...
    'FontName','Cambria')
    
%column 1
set(pq(1),'Position',[0.1 0.1 0.35 0.25],...
    'YTick',0:0.02:0.08,'YLim',[0 0.08],...
    'XTick',-2:1:2,'XLim',[-2 2])
set(pq(2),'Position',[0.1 0.42 0.35 0.25],...
    'YTick',0.82:0.02:0.9,'YLim',[0.8 0.9],...
    'XTick',0.6:0.1:1,'XLim',[0.6 1])
set(pq(3),'Position',[0.1 0.735 0.35 0.25],...
    'YTick',1.72:0.02:1.8,'YLim',[1.7 1.82],...
    'XTick',2:0.25:3,'XLim',[2 3])

%column 2
set(qp(1),'Position',[0.56 0.1 0.35 0.25],...
    'YTick',0:0.02:0.08,'YLim',[0 0.08],...
    'XTick',-0.4:0.2:0.4,'XLim',[-0.4 0.4])
set(qp(2),'Position',[0.56 0.42 0.35 0.25],...
    'YTick',0.82:0.02:0.9,'YLim',[0.8 0.9],...
    'XTick',-0.1:0.05:0.1,'XLim',[-0.1 0.1])
set(qp(3),'Position',[0.56 0.735 0.35 0.25],...
    'YTick',1.72:0.02:1.8,'YLim',[1.7 1.82],...
    'XTick',-0.1:0.005:-0.08,'XLim',[-0.1 -0.08])

export_fig(['c:\Users\bkn5\Projects\Mekong_W2015\Figures\Paper1\' 'VPmeanU&RS' name],'-jpg')

%%%%%%% plot #2: color plots
f2 = figure;
set(f2,'PaperOrientation','portrait',...
    'position',[400 100   1000   1000]);
set(gcf,'color','w','PaperPositionMode','auto')
colormap(redblue)
for i = 1:3
    ys = (dat.(fn{i}).rbcm)./ph;
    sbp1 = [5,3,1]; %subplot order
    pq(i) = subplot(3,2,sbp1(i));
    
    imagesc(dat.(fn{i}).T,ys,(dat.(fn{i}).U*100)')
    caxis([-2 6])
    set(gca,'Xlim',[dat.(fn{i}).T(1) dat.(fn{i}).T(end)],...
        'YDir','normal',...
        'FontSize',14,...
        'FontName','Cambria',...
        'LineWidth',1.5,...
        'TickLength',[0.02 0.02])
    hold off
    
    sbp2 = [6,4,2];
    qp(i) = subplot(3,2,sbp2(i));
    
    imagesc(dat.(fn{i}).T,ys,rss.(fn{i}).uw')
    caxis([-0.2 0.6])
    set(gca,'Xlim',[dat.(fn{i}).T(1) dat.(fn{i}).T(end)],...
        'YDir','normal',...
        'FontSize',14,...
        'FontName','Cambria',...
        'LineWidth',1.5,...
        'TickLength',[0.02 0.02])
    hold off
end
%plot adjustments
set([pq(2) pq(3) qp(2) qp(3)],'XTickLabel',[])
set(qp,'YTickLabel',[])
%column 1
set(pq(1),'Position',[0.1 0.23 0.36 0.22],'YTick',0:0.02:0.07,'YLim',[0 0.07])
set(pq(2),'Position',[0.1 0.485 0.36 0.22],'YTick',0.8:0.02:0.88,'YLim',[0.8 0.88])
set(pq(3),'Position',[0.1 0.75 0.36 0.22],'YTick',1.72:0.02:1.8,'YLim',[1.72 1.8])

%column 2
set(qp(1),'Position',[0.57 0.23 0.36 0.22],'YTick',0:0.02:0.07,'YLim',[0 0.07])
set(qp(2),'Position',[0.57 0.485 0.36 0.22],'YTick',0.8:0.02:0.88,'YLim',[0.8 0.88])
set(qp(3),'Position',[0.57 0.75 0.36 0.22],'YTick',1.72:0.02:1.8,'YLim',[1.72 1.8])

yl = ylabel(pq(2),'$$z/h_c$$');
set(yl,'Interpreter','latex','fontsize',14,'FontName','Cambria')

datetick(pq(1),'x','HH:MM:SS','keepticks','keeplimits')
datetick(qp(1),'x','HH:MM:SS','keepticks','keeplimits')

cb1 = colorbar('peer',pq(1),'southoutside');
cb2 = colorbar('peer',qp(1),'southoutside');

set(cb1,'position',[0.1 0.12 0.37 0.05],...
    'LineWidth',1.5,'FontSize',14,'FontName','Cambria')
x1 = xlabel(cb1,'U (cm/s)');

set(cb2,'position',[0.56 0.12 0.37 0.05],...
    'LineWidth',1.5,'FontSize',14,'FontName','Cambria')
x2 = xlabel(cb2,['$$\rho','\langle','\overline{u''v''}','   \rangle$$ (Pa)']);
set([x1 x2],'Interpreter','latex','fontsize',14,'FontName','Cambria')

export_fig(['c:\Users\bkn5\Projects\Mekong_W2015\Figures\Paper1\' 'VPmeanU&RS_colored' name],'-jpg')
