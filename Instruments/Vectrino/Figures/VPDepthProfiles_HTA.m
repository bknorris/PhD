%Calculate Depth profiles of mean velocities/Reyolds Stress for the
%Vectrino profilers, from HTA days 1-3

clear
datdir = 'c:\Users\bkn5\Projects\Mekong_W2015\DataAnalysis\Paper1\';
fname = dir([datdir 'HTAday*Vels.mat']);fname1 = {fname.name};
fname = dir([datdir '*day*RS.mat']);fname2 = {fname.name};
name = 'HTA';
ph = 49; %height of canopy in cm
heading = 10;
hab = [0.062 0.063 0.061;
    0.242 0.240 0.240;
    0.550 0.550 0.550]*100;
dn = {'day1';'day2';'day3'};
for i = 1:3
    load([datdir fname1{i}]) %velocity file
    load([datdir fname2{i}]) %reynolds stress file
    fn = fieldnames(HTA);
    if i == 1
        start = HTA.times.t1;stop = HTA.times.e1;
    elseif i == 2
        start = HTA.times.t2;stop = HTA.times.e2;
    elseif i == 3
        start = HTA.times.t3;stop = HTA.times.e3;
    end
    count = 1;
    for ii = 2:4
        %velocities
        idx = find(HTA.(fn{ii}).time >= start & HTA.(fn{ii}).time <= stop);
        u = HTA.(fn{ii}).y(idx,:); %along-shore
        v = HTA.(fn{ii}).x(idx,:); %cross-shore
        w = (HTA.(fn{ii}).z1(idx,:)+HTA.(fn{ii}).z2(idx,:))./2;
        t = HTA.(fn{ii}).time(idx,:);
        
        %rotate to the cross-shore direction
        rot = heading*pi/180;
        u = u.*(ones(size(u))*cos(rot)) + ...
            v.*(ones(size(v))*sin(rot));
        v = -u.*(ones(size(u))*sin(rot)) + ...
            v.*(ones(size(v))*cos(rot));
        
        h = hab(i,count); %inst hab (cm)
        rb = HTA.(fn{ii}).rb;dat.(dn{i}).(fn{ii}).rbcm = h-rb(1:end)*100;
        sr = 50;
        intv = 10; %minutes
        avt = intv*sr*60;
        ind = [1 avt:avt:length(u)];
        
        for j = 1:length(ind)-1
            Uz = u(ind(j):ind(j+1),:);
            dat.(dn{i}).(fn{ii}).U(j,:) = nanmean(Uz);
            dat.(dn{i}).(fn{ii}).T(j,:) = t(ind(j));
        end
        dat.(dn{i}).(fn{ii}).uw = rss.(fn{ii}).uw;
        count = count+1;
    end
end

%%%%%%%% Plotting Routine
fn = fieldnames(dat.day1);
f1 = figure;
set(f1,'PaperOrientation','portrait',...
    'position',[400 100   800   1000]);
set(gcf,'color','w','PaperPositionMode','auto')

for i = 1:3
    c = [0.6 0.6 0.6;0.4 0.4 0.4;0 0 0];
    symb = {'s';'^';'o'};
    sbp1 = [5,3,1]; %subplot order
    pq(i) = subplot(3,2,sbp1(i));
    for ii = 1:3
        mVel = mean(dat.(dn{i}).(fn{ii}).U);Velstd = std(dat.(dn{i}).(fn{ii}).U);
        ys = (dat.(dn{i}).(fn{ii}).rbcm)./ph;
        
        zby = linspace(0,2,length(ys));
        zbx = zeros(length(zby),1);
        plot(zbx,zby,'--k','LineWidth',1.5), hold on
        %         h1 = herrorbar(mVel(1:3:end)*100,ys(1:3:end),2*100*Velstd(1:3:end));
        %         set(h1,'LineWidth',1.5,'Color','k')
        plot(mVel(1:3:end)*100,ys(1:3:end),'Color',c(ii,:),'LineWidth',2,...
            'marker',symb{ii},'markersize',8,'markerfacecolor',[1 1 1])
    end
    xl = xlabel('U (cm/s)');
    set(xl,'Interpreter','latex','fontsize',14,'FontName','Cambria')
    yl = ylabel('$$z/h_c$$');
    set(yl,'Interpreter','latex','fontsize',14,'FontName','Cambria')
    grid on
    set(gca,'GridLineStyle',':')
    if sbp1(i) == 1
        xs = linspace(-5,10,10);
        ys = ones(10,1);
        plot(xs,ys,'--k','LineWidth',1.5)
    end
end

for i = 1:3
    sbp2 = [6,4,2]; %subplot order
    qp(i) = subplot(3,2,sbp2(i));
    for ii = 1:3
        mRS = mean(dat.(dn{i}).(fn{ii}).uw);RSstd = std(dat.(dn{i}).(fn{ii}).uw);
        ys = (dat.(dn{i}).(fn{ii}).rbcm)./ph;
        plot(zbx,zby,'--k','LineWidth',1.5), hold on
        %         h2 = herrorbar(mRS(1:3:end),ys(1:3:end),2*RSstd(1:3:end));
        %         set(h2,'LineWidth',1.5,'Color','k')
        pp(ii) = plot(mRS(1:3:end),ys(1:3:end),'Color',c(ii,:),'LineWidth',2,...
            'marker',symb{ii},'markersize',8,'markerfacecolor',[1 1 1]);
    end
    xl = xlabel(['$$\rho','\langle','\overline{u''v''}','    \rangle$$ (Pa)']);
    set(xl,'Interpreter','latex','fontsize',14,'FontName','Cambria')
    yl = ylabel('$$z/h_c$$');
    set(yl,'Interpreter','latex','fontsize',14,'FontName','Cambria')
    grid on
    set(gca,'GridLineStyle',':')
    if sbp2(i) == 2
        xs = linspace(-0.5,0.5,10);
        ys = ones(10,1);
        plot(xs,ys,'--k','LineWidth',1.5)
    end
end
%legend
leg = legend(pp,'x = -10cm','x = 10cm','x = 20cm');
set(leg,'position',[0.82 0.58 0.05 0.05],'LineWidth',1.5)
%plot adjustments
set([pq qp],'LineWidth',1.5,'FontSize',14,...
    'FontName','Cambria')
    
%column 1
set(pq(1),'Position',[0.1 0.1 0.35 0.25],...
    'YTick',0:0.02:0.06,'YLim',[0 0.06],...
    'XTick',-2:1:2,'XLim',[-2 2])
set(pq(2),'Position',[0.1 0.42 0.35 0.25],...
    'YTick',0.36:0.02:0.42,'YLim',[0.34 0.42],...
    'XTick',-0.8:0.4:0.8,'XLim',[-0.8 0.8])
set(pq(3),'Position',[0.1 0.735 0.35 0.25],...
    'YTick',0.98:0.02:1.04,'YLim',[0.96 1.05],...
    'XTick',3.6:0.2:4.4,'XLim',[3.6 4.4])

%column 2
set(qp(1),'Position',[0.56 0.1 0.35 0.25],...
    'YTick',0:0.02:0.06,'YLim',[0 0.06],...
    'XTick',-0.4:0.2:0.4,'XLim',[-0.4 0.4])
set(qp(2),'Position',[0.56 0.42 0.35 0.25],...
    'YTick',0.36:0.02:0.4,'YLim',[0.34 0.42],...
    'XTick',-0.2:0.1:0.2,'XLim',[-0.2 0.2])
set(qp(3),'Position',[0.56 0.735 0.35 0.25],...
    'YTick',0.98:0.02:1.04,'YLim',[0.96 1.05],...
    'XTick',-0.4:0.1:0,'XLim',[-0.4 0])

export_fig(['c:\Users\bkn5\Projects\Mekong_W2015\Figures\Paper1\' 'VPmeanU&RS' name],'-jpg')

