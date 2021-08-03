%TKE dissipation rate time series for various deployments, this will be for the 2014 experiments to start
clear

tkedir = 'd:\Projects\Mekong_F2014\DataAnalysis\Paper2\TKE\';
tkefiles = dir(tkedir);tkefiles = {tkefiles.name};
%Julia wants to see F2F (all instruments) and MA1 (which has the highest
%TKE dissipation rate from the VegDensityTKE figure). Let's create a time
%series plot of time/E, color markers by significant wave height.
%average TKE estimates (beam1+beam3)/2


%flats 2 forest figure:
load([tkedir tkefiles{5}])
load('d:\Projects\Mekong_F2014\DataAnalysis\Paper2\Spectra\WaveStatsHR3_F2F_1.mat')
start = Stat.vpro1.time(1);stop = Stat.vpro1.time(end); %define start and end times to crop wvstats
ind = find(wvstats.time >= start & wvstats.time <= stop);ind(9) = 37; %only for F2F_2
Hsig = sqrt(2)*wvstats.hrmsuv(ind);

%crop wave statistics file to the duration of the TKE estimates
vp1 = (Stat.vpro1.beam1.E+Stat.vpro1.beam3.E)./2;
vp2 = (Stat.vpro2.beam1.E+Stat.vpro2.beam3.E)./2;
vp3 = (Stat.vpro3.beam1.E+Stat.vpro3.beam3.E)./2;
[~,m] = size(vp1);Hsig = Hsig(1:m);
time = (1:1:m)*10;
bins = [5 10 15 20 25];
cmap = flipud([0,242,255;0,193,250;0,149,237;0,97,219;0,0,128]./256);
marker = {'o';'s';'^';'d';'>'};


f1 = figure(1);
set(f1,'PaperOrientation','portrait',...
    'position',[400 200   1400   500]);
set(gcf,'color','w','PaperPositionMode','auto')
%color by Hsig
crange = 0:0.05:0.4;
wvmap = linspecer(length(crange),'sequential');
sp(1) = subplot(131);
for i = 1:length(bins)
    plot(time,log10(vp1(bins(i),:)),'LineWidth',3,'Color',cmap(i,:)), hold on
    for ii = 1:m
        tmp = abs(crange-Hsig(ii));
        [~,idx] = min(tmp);
        pp(i)= plot(time(ii),log10(vp1(bins(i),ii)),'Marker',marker{i},'LineWidth',1.5,...
            'MarkerFaceColor',wvmap(idx,:),'MarkerEdgeColor','k','MarkerSize',10);
    end
end
hold off
sp(2) = subplot(132);
for i = 1:length(bins)
    plot(time,log10(vp2(bins(i),:)),'LineWidth',3,'Color',cmap(i,:)), hold on
    for ii = 1:m
        tmp = abs(crange-Hsig(ii));
        [~,idx] = min(tmp);
        plot(time(ii),log10(vp2(bins(i),ii)),'Marker',marker{i},'LineWidth',1.5,...
            'MarkerFaceColor',wvmap(idx,:),'MarkerEdgeColor','k','MarkerSize',10)
    end
end
hold off
sp(3) = subplot(133);
for i = 1:length(bins)
    plot(time,log10(vp3(bins(i),:)),'LineWidth',3,'Color',cmap(i,:)), hold on
    for ii = 1:m
        tmp = abs(crange-Hsig(ii));
        [~,idx] = min(tmp);
        plot(time(ii),log10(vp3(bins(i),ii)),'Marker',marker{i},'LineWidth',1.5,...
            'MarkerFaceColor',wvmap(idx,:),'MarkerEdgeColor','k','MarkerSize',10)
    end
end
hold off
%colorbar
colormap(wvmap);
cb = colorbar;
caxis([0 0.4])
%global plot adjustments
set(sp,'Xlim',[0 time(end)],'XTick',0:20:100,'Ylim',[-6.5 -2.5],...
    'LineWidth',1.5,'FontSize',14,'FontName','Cambria',...
    'Box','on')
set([sp(2) sp(3)],'YTickLabel',[])
set(cb,'LineWidth',1.5,'FontSize',14,'FontName','Cambria')
%positioning
set(sp(1),'position',[0.08 0.15 0.25 0.75])
set(sp(2),'position',[0.365 0.15 0.25 0.75])
set(sp(3),'position',[0.65 0.15 0.25 0.75])
set(cb,'position',[0.92 0.15 0.01 0.75])
%labeling
ylabel(sp(1),'Log_1_0(\epsilon) (W/kg)')
xlabel(sp(2),'Time after submergence (min)')
ylabel(cb,'H_s (m)')
%legend
leg = legend(pp,{'Bin 5';'Bin 10';'Bin 15';'Bin 20';'Bin 25'});
set(leg,'position',[0.5 0.25 0.1 0.1],'box','off')
% figdir = 'c:\Users\bkn5\Projects\Mekong_F2014\Figures\Paper2\';
% export_fig([figdir 'F2F_300914_ts'],'-jpeg','-nocrop')

%MA1 figure:
load([tkedir tkefiles{10}])
load('c:\Users\bkn5\Projects\Mekong_F2014\DataAnalysis\Paper2\Spectra\WaveStatsHR3_MA1.mat')
start = Stat.vpro1.time(1);stop = Stat.vpro1.time(end); %define start and end times to crop wvstats
ind = find(wvstats.time >= start & wvstats.time <= stop);ind(6) = 8;
Hsig = sqrt(2)*wvstats.hrmsuv(ind);

%crop wave statistics file to the duration of the TKE estimates
vp1 = (Stat.vpro1.beam1.E+Stat.vpro1.beam3.E)./2;
vp2 = (Stat.vpro2.beam1.E+Stat.vpro2.beam3.E)./2;
vp3 = (Stat.vpro3.beam1.E+Stat.vpro3.beam3.E)./2;
[~,m] = size(vp1);Hsig = Hsig(1:m);
time = (1:1:m)*10;
bins = [5 10 15 20 25];
cmap = flipud([0,242,255;0,193,250;0,149,237;0,97,219;0,0,128]./256);
marker = {'o';'s';'^';'d';'>'};


f2 = figure(2);
set(f2,'PaperOrientation','portrait',...
    'position',[400 200   1400 500]);
set(gcf,'color','w','PaperPositionMode','auto')
%color by Hsig
crange = 0:0.05:0.3;
wvmap = linspecer(length(crange),'sequential');
sp(1) = subplot(131);
for i = 1:length(bins)
    plot(time,log10(vp1(bins(i),:)),'LineWidth',3,'Color',cmap(i,:)), hold on
    for ii = 1:m
        tmp = abs(crange-Hsig(ii));
        [~,idx] = min(tmp);
        pp(i)= plot(time(ii),log10(vp1(bins(i),ii)),'Marker',marker{i},'LineWidth',1.5,...
            'MarkerFaceColor',wvmap(idx,:),'MarkerEdgeColor','k','MarkerSize',10);
    end
end
hold off
sp(2) = subplot(132);
for i = 1:length(bins)
    plot(time,log10(vp2(bins(i),:)),'LineWidth',3,'Color',cmap(i,:)), hold on
    for ii = 1:m
        tmp = abs(crange-Hsig(ii));
        [~,idx] = min(tmp);
        plot(time(ii),log10(vp2(bins(i),ii)),'Marker',marker{i},'LineWidth',1.5,...
            'MarkerFaceColor',wvmap(idx,:),'MarkerEdgeColor','k','MarkerSize',10)
    end
end
hold off
sp(3) = subplot(133);
for i = 1:length(bins)
    plot(time,log10(vp3(bins(i),1:m)),'LineWidth',3,'Color',cmap(i,:)), hold on
    for ii = 1:m
        tmp = abs(crange-Hsig(ii));
        [~,idx] = min(tmp);
        plot(time(ii),log10(vp3(bins(i),ii)),'Marker',marker{i},'LineWidth',1.5,...
            'MarkerFaceColor',wvmap(idx,:),'MarkerEdgeColor','k','MarkerSize',10)
    end
end
hold off
%colorbar
colormap(wvmap);
cb = colorbar;
caxis([0 0.3])
%global plot adjustments
set(sp,'Xlim',[0 time(end)],'XTick',0:20:100,'Ylim',[-6.5 -2.5],...
    'LineWidth',1.5,'FontSize',14,'FontName','Cambria',...
    'Box','on')
set([sp(2) sp(3)],'YTickLabel',[])
set(cb,'LineWidth',1.5,'FontSize',14,'FontName','Cambria')
%positioning
set(sp(1),'position',[0.08 0.15 0.25 0.75])
set(sp(2),'position',[0.365 0.15 0.25 0.75])
set(sp(3),'position',[0.65 0.15 0.25 0.75])
set(cb,'position',[0.92 0.15 0.01 0.75])
%labeling
ylabel(sp(1),'Log_1_0(\epsilon) (W/kg)')
xlabel(sp(2),'Time after submergence (min)')
ylabel(cb,'H_s (m)')
%legend
leg = legend(pp,{'Bin 5';'Bin 10';'Bin 15';'Bin 20';'Bin 25'});
set(leg,'position',[0.5 0.25 0.1 0.1],'box','off')
figdir = 'c:\Users\bkn5\Projects\Mekong_F2014\Figures\Paper2\';
export_fig([figdir 'FSS_1_250914_ts'],'-jpeg','-nocrop')