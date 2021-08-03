%structure function statistics figures, version 2:
%1. Timeseries of b (Struct. fun slope) and time for a single depth bin, say
%   bin 15. Color by r^2.
%2. Timeseries of the error in epsilon, for all instruments over the three
%   days. Color by waveheight. 

clear
savedatdir = 'c:\Users\bkn5\Projects\Mekong_W2015\DataAnalysis\Turbulence\Paper1\';
figdir = 'c:\Users\bkn5\Projects\Mekong_W2015\Figures\Turbulence\Statistics\';
load('C:\Users\bkn5\Projects\Mekong_W2015\DataAnalysis\Spectra\Paper1\WaveStatsV5108_HTA.mat') %wave statistics
load('C:\Users\bkn5\Projects\Mekong_W2015\DataAnalysis\Paper1\HTAtimes.mat')
fname = dir([savedatdir 'HTA*']);fname = {fname.name};
dn = {'day1';'day2';'day3'};
%%%%% Data conglomerator
for i = 1:3                                                                 %loop through days
    load([savedatdir fname{i}])
    fn = fieldnames(Stat);
    for ii = 1:3                                                            %loop through instruments
        %calculate average fits between beams 1 and 3
        bfit = (Stat.(fn{ii}).beam1.bfit+Stat.(fn{ii}).beam3.bfit)./2;
        rsq = (Stat.(fn{ii}).beam1.rsq+Stat.(fn{ii}).beam3.rsq)./2;
        pval = (Stat.(fn{ii}).beam1.pval+Stat.(fn{ii}).beam3.pval)./2;
        
        %get error in epsilon
        epsilon = (Stat.(fn{ii}).beam1.E+Stat.(fn{ii}).beam3.E)./2;
        epserr = (Stat.(fn{ii}).beam1.Emaxerr+Stat.(fn{ii}).beam3.Emaxerr)./2;
        
        %no times in Stat file. Use times from the HTAtimes file. Interval
        %step in Stat file is 10 minutes.
        if i == 1
            start = HTA.times.t1;stop = HTA.times.e1;
        elseif i == 2
            start = HTA.times.t2;stop = HTA.times.e2;
        elseif i == 3
            start = HTA.times.t3;stop = HTA.times.e3;
        end
        [~,m] = size(Stat.(fn{ii}).beam1.E);
        intv = datenum(0,0,0,0,10,0);
        time = start:intv:stop;
        idx = zeros(m,1);
        for iii = 1:m
            tmp = abs(wvstats.time-time(iii));
            [~,idx(iii)] = min(tmp);
        end
        
        %build data structure
        dat.(dn{i}).(fn{ii}).Time = time(1:m);
        dat.(dn{i}).(fn{ii}).bfit = bfit;
        dat.(dn{i}).(fn{ii}).pval = pval;
        dat.(dn{i}).(fn{ii}).rsq = rsq;
        dat.(dn{i}).(fn{ii}).E = epsilon;
        dat.(dn{i}).(fn{ii}).Err = epserr;
        dat.(dn{i}).(fn{ii}).hrms = wvstats.hrmsp(idx);
    end
end

%%%%% Plot Routine
f1 = figure(1);
set(f1,'PaperOrientation','portrait',...
    'position',[400 100   800   1000]);
set(gcf,'color','w','PaperPositionMode','auto')
gap = datenum(0,0,0,0,5,0);     %adjusts the xaxis lims so symbols don't overlap the plot box
%first column will be the b timeseries, colored by r^2
rc = 0.95:0.001:0.999;                                                %possible values for r^2
c = linspecer(length(rc),'sequential');colormap(c);
c2 = [0.6 0.6 0.6;0.4 0.4 0.4;0 0 0];                               %line plot colorspec
for i = 1:3                                                         %days
    symb = {'s';'^';'o'};
    sbp1 = [5,3,1];                                                 %subplot order
    pq(i) = subplot(3,2,sbp1(i));
    for ii = 1:3                                                    %instruments
        [~,m] = size(dat.(dn{i}).(fn{ii}).Time);
        if i == 1 && ii == 3    %day 1, VP3 bin 15 is below the bed surface, use bin 7 consistent with other analyses
            bin = 7;
        else
            bin = 15;
        end
        plot(dat.(dn{i}).(fn{ii}).Time,dat.(dn{i}).(fn{ii}).bfit(bin,:),...
            'LineWidth',1.5,'Color',c2(ii,:)), hold on
        for iii = 1:m
            tmp = abs(rc-dat.(dn{i}).(fn{ii}).rsq(bin,iii));
            [~,idx] = min(tmp);
            plot(dat.(dn{i}).(fn{ii}).Time(iii),dat.(dn{i}).(fn{ii}).bfit(bin,iii),...
            'Marker',symb{ii},'MarkerSize',8,'MarkerFaceColor',c(idx,:),...
            'LineWidth',1.5,'MarkerEdgeColor',c2(ii,:))
        end
        
    end
    grid on
    set(gca,'GridLineStyle',':')
    xs = linspace(dat.(dn{i}).vpro1.Time(1)-gap,dat.(dn{i}).vpro1.Time(end)+gap,10);
    ys = (2/3)*ones(10,1);
    plot(xs,ys,'--k','LineWidth',1.5)
    yl = ylabel('r^b slope fits');
    set(yl,'fontsize',12,'FontName','Cambria')
end

%second column will be the error timeseries, colored by waveheight    
wh = min(dat.day1.vpro1.hrms):0.01:max(dat.day3.vpro1.hrms); %possible values for waveheight
c = linspecer(length(wh),'sequential');colormap(c);
for i = 1:3                                                         %days
    sbp2 = [6,4,2];                                                 %subplot order
    qp(i) = subplot(3,2,sbp2(i));
    for ii = 1:3                                                    %instruments
        [~,m] = size(dat.(dn{i}).(fn{ii}).Time);
        if i == 1 && ii == 3    %day 1, VP3 bin 15 is below the bed surface, use bin 7 consistent with other analyses
            bin = 7;
        else
            bin = 15;
        end
        norm = range(dat.(dn{i}).(fn{ii}).hrms);                    %normalize error estimate by the data range
        plot(dat.(dn{i}).(fn{ii}).Time,dat.(dn{i}).(fn{ii}).Err(bin,:)./norm,...
            'LineWidth',1.5,'Color',c2(ii,:)), hold on
        for iii = 1:m
            tmp = abs(wh-dat.(dn{i}).(fn{ii}).hrms(iii));
            [~,idx] = min(tmp);
            plot(dat.(dn{i}).(fn{ii}).Time(iii),dat.(dn{i}).(fn{ii}).Err(bin,iii)/norm,...
            'Marker',symb{ii},'MarkerSize',8,'MarkerFaceColor',c(idx,:),...
            'LineWidth',1.5,'MarkerEdgeColor',c2(ii,:)),hold on
        end
    end
    yl = ylabel('$$\hat{\epsilon}_{err}$$');
    set(yl,'Interpreter','latex','fontsize',16,'FontName','Cambria')
    grid on
    set(gca,'GridLineStyle',':')
end
hold off

%labeling
xl = xlabel(qp(1),'Time on 07/03/2015');
set(xl,'fontsize',12,'FontName','Cambria')
xl = xlabel(pq(1),'Time on 07/03/2015');
set(xl,'fontsize',12,'FontName','Cambria')
xl = xlabel(qp(2),'Time on 08/03/2015');
set(xl,'fontsize',12,'FontName','Cambria')
xl = xlabel(pq(2),'Time on 08/03/2015');
set(xl,'fontsize',12,'FontName','Cambria')
xl = xlabel(qp(3),'Time on 10/03/2015');
set(xl,'fontsize',12,'FontName','Cambria')
xl = xlabel(pq(3),'Time on 10/03/2015');
set(xl,'fontsize',12,'FontName','Cambria')

%colorbars
cb1 = colorbar('peer',pq(1),'southoutside');caxis(pq(1),[0.95 0.99]) 
set(cb1,'position',[0.1 0.08 0.37 0.02],...
    'LineWidth',1.5,'FontSize',12,'FontName','Cambria')
xlabel(cb1,'r^2');

cb2 = colorbar('peer',qp(1),'southoutside');caxis(qp(1),[0 0.6])
set(cb2,'position',[0.6 0.08 0.37 0.02],...
    'LineWidth',1.5,'FontSize',12,'FontName','Cambria',...
    'XTick',0:0.2:0.6)
xlabel(cb2,'H_r_m_s (m)');

%global plot adjustments
set([pq qp],'LineWidth',1.5,'FontSize',12,'FontName','Cambria')
step = [4 6 2];
for i = 1:3
    t1 = dat.(dn{i}).vpro1.Time(1);t2 = dat.(dn{i}).vpro1.Time(end);
    set(pq(i),'Xlim',[t1-gap t2+gap],'XTick',t1:intv*step(i):t2,...
        'Ylim',[0 1.6],'YTick',0:0.5:1.5)
    set(qp(i),'Xlim',[t1-gap t2+gap],'XTick',t1:intv*step(i):t2,...
        'Ylim',[0 1E-2],'YTick',0:0.0025:1E-2)
    datetick(pq(i),'x','HH:MM','keepticks','keeplimits')
    datetick(qp(i),'x','HH:MM','keepticks','keeplimits')
end
set([qp(2) qp(3)],'YLim',[0 3E-3],'YTick',0:0.00075:3E-3)

%positioning
%column 1
set(pq(1),'Position',[0.1 0.18 0.37 0.22])
set(pq(2),'Position',[0.1 0.47 0.37 0.22])
set(pq(3),'Position',[0.1 0.75 0.37 0.22])

%column 2
set(qp(1),'Position',[0.6 0.18 0.37 0.22])
set(qp(2),'Position',[0.6 0.47 0.37 0.22])
set(qp(3),'Position',[0.6 0.75 0.37 0.22])

export_fig(['c:\Users\bkn5\Projects\Mekong_W2015\Figures\Paper1\' 'VP_SFts_HTA'],'-jpg','-nocrop','-m1')

%plot pvalue of the timeseries for each day
f2 = figure(2);
set(f2,'PaperOrientation','portrait',...
    'position',[400 100   1000   400]);
set(gcf,'color','w','PaperPositionMode','auto')
wh = min(dat.day1.vpro1.hrms):0.01:max(dat.day3.vpro1.hrms); %possible values for waveheight
c = linspecer(length(wh),'sequential');colormap(c);
c2 = [0.6 0.6 0.6;0.4 0.4 0.4;0 0 0];                               %line plot colorspec
gap = datenum(0,0,0,0,5,0);     %adjusts the xaxis lims so symbols don't overlap the plot box
for i = 1:3                                                         %days
    symb = {'s';'^';'o'};
    qp(i) = subplot(1,3,i);                                         %subplots
    for ii = 1:3                                                    %instruments
        [~,m] = size(dat.(dn{i}).(fn{ii}).Time);
        if i == 1
        bin = 1;
        else
            bin = 13;
        end
        plot(dat.(dn{i}).(fn{ii}).Time,dat.(dn{i}).(fn{ii}).pval(bin,:),...
            'LineWidth',1.5,'Color',c2(ii,:)), hold on
        for iii = 1:m
            tmp = abs(wh-dat.(dn{i}).(fn{ii}).hrms(iii));
            [~,idx] = min(tmp);
            plot(dat.(dn{i}).(fn{ii}).Time(iii),dat.(dn{i}).(fn{ii}).pval(bin,iii),...
            'Marker',symb{ii},'MarkerSize',8,'MarkerFaceColor',c(idx,:),...
            'LineWidth',1.5,'MarkerEdgeColor',c2(ii,:)),hold on
        end
    end
    ylabel('P-value','fontsize',12,'FontName','Cambria')
    grid on
    set(gca,'GridLineStyle',':')
end
hold off

%labeing
xl = xlabel(qp(1),'Time on 07/03/2015');
set(xl,'fontsize',12,'FontName','Cambria')
xl = xlabel(qp(2),'Time on 08/03/2015');
set(xl,'fontsize',12,'FontName','Cambria')
xl = xlabel(qp(3),'Time on 10/03/2015');
set(xl,'fontsize',12,'FontName','Cambria')

%colorbar
cb2 = colorbar('peer',qp(3),'eastoutside');caxis(qp(3),[0 0.6])
set(cb2,'position',[0.91 0.15 0.02 0.8],...
    'LineWidth',1.5,'FontSize',12,'FontName','Cambria',...
    'YTick',0:0.1:0.6)
ylabel(cb2,'H_r_m_s (m)');

%global plot adjustments
set(qp,'LineWidth',1.5,'FontSize',10,'FontName','Cambria')
step = [4 6 2];
for i = 1:3
    t1 = dat.(dn{i}).vpro1.Time(1);t2 = dat.(dn{i}).vpro1.Time(end);
    set(qp(i),'Xlim',[t1-gap t2+gap],'XTick',t1:intv*step(i):t2,...
        'Ylim',[0 0.05],'YTick',0:0.01:0.05)
    datetick(qp(i),'x','HH:MM','keepticks','keeplimits')
end
set([qp(2) qp(3)],'YLim',[0 0.05],'YTick',0:0.01:0.05)

%positioning
set(qp(1),'Position',[0.07 0.15 0.22 0.8])
set(qp(2),'Position',[0.37 0.15 0.22 0.8])
set(qp(3),'Position',[0.67 0.15 0.22 0.8])

export_fig(['c:\Users\bkn5\Projects\Mekong_W2015\Figures\Paper1\' 'HTA_pvals'],'-jpg','-nocrop','-m1')
