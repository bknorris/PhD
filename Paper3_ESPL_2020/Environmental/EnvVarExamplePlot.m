%Plot overview figure of bdmed, depth, Hs, umag for a single experiment
%(flood, ebb);
clear,close all
load('e:\Mekong_W2015\DataAnalysis\Paper3\CmbData\11-03-15\AllEventData.mat');
dat.flood = data;clear data
load('e:\Mekong_W2015\DataAnalysis\Paper3\CmbData\12-03-15\AllEventData.mat');
dat.ebb = data;clear data
fn = fieldnames(dat);ffn = fieldnames(dat.flood);
sfdir = 'd:\GradSchool\DataAnalysis\Paper3\Figures\';
flds = {'Mudflat';'Fringe';'Forest'};
f1 = figure(1);
set(f1,'PaperOrientation','portrait',...
    'position',[500 10   900   700],...
    'renderer','painters');
vars = {'depth';'Hs';'umag';'bdist'};
w = [1 3 5 7;2 4 6 8];
cl = [0.7 0.7 0.7;0.4 0.4 0.4;0.1 0.1 0.1];
sp = zeros(2,4);
pp = zeros(3,1);
for i = 1:2
    for j = 1:length(vars)
        sp(w(i,j)) = subplot(4,2,w(i,j));
        hold on
        for k = 1:3
            time = dat.(fn{i}).(ffn{k}).time;
            y1 = dat.(fn{i}).(ffn{k}).(vars{j});
            if strcmp(vars{j},'umag')
                y1 = fastsmooth(y1,50,3,1);
            end
            if strcmp(vars{j},'bdist')
                    y1 = y1-0.061;
                if k == 1
                plot(time,zeros(length(time),1),'-k',...
                    'linewidth',1)
                end
            end
            pp(k) = plot(time,y1,'linewidth',1.5,...
                'color',cl(k,:));
        end
    end
end
leg = legend(pp,flds);
set(leg,'position',[0.9 0.89 0.05 0.05])
t1 = dat.flood.vpro1.time;t2 = dat.ebb.vpro1.time;
tstep = datenum(0,0,0,0,30,0);
set([sp(1),sp(2),sp(3),sp(4),sp(5),sp(6)],'xticklabel',[])
set([sp(2),sp(4),sp(6),sp(8)],'yticklabel',[])
set(sp(1),'xlim',[t1(1) t1(end)],'xtick',t1(1):tstep:t1(end),'ylim',[0 1.5],'position',[0.12 0.78 0.4 0.18])
set(sp(2),'xlim',[t2(1) t2(end)],'xtick',t2(1):tstep:t2(end),'ylim',[0 1.5],'position',[0.55 0.78 0.4 0.18])
set(sp(3),'xlim',[t1(1) t1(end)],'xtick',t1(1):tstep:t1(end),'ylim',[0 0.6],'position',[0.12 0.56 0.4 0.18])
set(sp(4),'xlim',[t2(1) t2(end)],'xtick',t2(1):tstep:t2(end),'ylim',[0 0.6],'position',[0.55 0.56 0.4 0.18])
set(sp(5),'xlim',[t1(1) t1(end)],'xtick',t1(1):tstep:t1(end),'ylim',[0 0.3],'position',[0.12 0.335 0.4 0.18])
set(sp(6),'xlim',[t2(1) t2(end)],'xtick',t2(1):tstep:t2(end),'ylim',[0 0.3],'position',[0.55 0.335 0.4 0.18])
set(sp(7),'xlim',[t1(1) t1(end)],'xtick',t1(1):tstep:t1(end),'ylim',[-0.04 0.04],...
    'ytick',-0.04:0.02:0.04,'position',[0.12 0.12 0.4 0.18])
set(sp(8),'xlim',[t2(1) t2(end)],'xtick',t2(1):tstep:t2(end),'ylim',[-0.04 0.04],...
    'ytick',-0.04:0.02:0.04,'position',[0.55 0.12 0.4 0.18])
datetick(sp(7),'x','HH:MM','keepticks','keeplimits')
datetick(sp(8),'x','HH:MM','keepticks','keeplimits')
xlabel(sp(7),['Time on ' datestr(t1(1),'dd-mm-yy')])
xlabel(sp(8),['Time on ' datestr(t2(1),'dd-mm-yy')])
title(sp(1),'Flood Tide')
title(sp(2),'Ebb Tide')
ylabel(sp(1),'Water Depth [m]')
ylabel(sp(3),'H_s [m]')
ylabel(sp(5),'$\overline{u}\quad [ms^{-1}]$','interpreter','latex')
ylabel(sp(7),'BLE [m]')
prettyfigures('text',12,'labels',13,'box',1)
export_fig([sfdir 'EnvVarExamp_fld_ebb'],'-pdf','-nocrop')