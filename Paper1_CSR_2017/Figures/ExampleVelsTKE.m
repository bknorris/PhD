%plot average velocities, TKE dissipation rate for F2F2_2 as an example of
%the dataset for Paper 1.

clear
close all

veldir = 'd:\Projects\Mekong_W2015\DataAnalysis\Paper1\AveragedVelocities\Original\';
tkedir = 'd:\Projects\Mekong_W2015\DataAnalysis\Paper1\TKE\Vertical\';
load([veldir 'F2F2_2Vave_bdadj.mat'])
load([tkedir 'F2F2_2TKE_bdadj.mat'])

start = datenum(2015,03,06,13,40,00);
stop = datenum(2015,03,06,15,40,00);
twm = datenum(0,0,0,0,20,0);ten = datenum(0,0,0,0,10,0);
for k = 1:3
        %%%Average TKE to the same length as Avgs
        fn = fieldnames(Avgs);
%         n = length(Avgs.(fn{k}).time);
%         beam1 = zeros(30,n);beam2 = zeros(30,n);beam3 = zeros(30,n);beam4 = zeros(30,n);
%         tt = zeros(1,n);
%         for kk = 1:n-1
%             td = find(Stat.(fn{k}).time >= Avgs.(fn{k}).time(kk) & Stat.(fn{k}).time <= Avgs.(fn{k}).time(kk+1));
%             beam1(1:30,kk) = nanmean(Stat.(fn{k}).beam1.E(:,td),2);
%             beam2(1:30,kk) = nanmean(Stat.(fn{k}).beam2.E(:,td),2);
%             beam3(1:30,kk) = nanmean(Stat.(fn{k}).beam3.E(:,td),2);
%             beam4(1:30,kk) = nanmean(Stat.(fn{k}).beam4.E(:,td),2);
%             tt(kk) = Stat.(fn{k}).time(td(1));
%         end
%         tt(end) = Stat.(fn{k}).time(td(end));
        tt = Stat.(fn{k}).time;
%         tke = (beam1+beam2+beam3+beam4)./4;
        z1 = Stat.(fn{k}).z1.E;z2 = Stat.(fn{k}).z2.E;
        tke = (z1+z2)./2;
        ind = find(tt >= start & tt <= stop);
        tt = tt(ind);tke = tke(:,ind);
        Stat.(fn{k}).tke = tke;Stat.(fn{k}).tt = tt;
        ind = find(Avgs.(fn{k}).time >= start & Avgs.(fn{k}).time <= stop);
        Avgs.(fn{k}).av = Avgs.(fn{k}).x(:,ind);
        Avgs.(fn{k}).tt = Avgs.(fn{k}).time(ind);
end


%Plot figures
f1 = figure(1);
set(f1,'PaperOrientation','portrait',...
    'position',[400 100   1100 900]);
set(gcf,'color','w','PaperPositionMode','auto')

%%%Average velocities
c = cbrewer('seq','Greys',100,'pchip');
c(2:101,:) = c(1:100,:);
c(1,:) = [0.9 0.9 0.9]; %set -inf value color
colormap(c)
z = linspace(-10,20,30);z = fliplr(z);

sp(1) = subplot(331);
tt = Avgs.vpro1.tt;av = Avgs.vpro1.av;
av(av == 0) = -inf;
imagesc(tt,z,av)
caxis([-0.04 0.04])
set(gca,'ydir','normal')
freezeColors

sp(2) = subplot(332);
tt = Avgs.vpro2.tt;av = Avgs.vpro2.av;
av(av == 0) = -inf;
imagesc(tt,z,av)
caxis([-0.04 0.04])
set(gca,'ydir','normal')
freezeColors

sp(3) = subplot(333);
tt = Avgs.vpro3.tt;av = Avgs.vpro3.av;
av(av == 0) = -inf;
imagesc(tt,z,av)
caxis([-0.04 0.04])
set(gca,'ydir','normal')

cb1 = colorbar;
freezeColors
cdh = get(cb1,'children');
map = get(cdh,'cdata');
set(cdh,'cdata',map(2:end));
cb1 = cbfreeze(cb1);


%%%Epsilon estimates
c = cbrewer('seq','Greys',100,'pchip');
c(2:101,:) = c(1:100,:);
c(1,:) = [0.9 0.9 0.9]; %set -inf value color
colormap(c)

sp(4) = subplot(334);
tt = Stat.vpro1.tt;eps = Stat.vpro1.tke;
eps(eps == 0) = -inf;
eps(isnan(eps)) = -inf;
imagesc(tt,z,eps)
caxis([1E-5 3E-4])
freezeColors
datetickzoom('x','HH:MM','keepticks','keeplimits')
set(gca,'ydir','normal')

sp(5) = subplot(335);
tt = Stat.vpro2.tt;eps = Stat.vpro2.tke;
eps(eps == 0) = -inf;
eps(isnan(eps)) = -inf;
imagesc(tt,z,eps)
caxis([1E-5 3E-4])
freezeColors
datetickzoom('x','HH:MM','keepticks','keeplimits')
set(gca,'ydir','normal')

sp(6) = subplot(336);
tt = Stat.vpro3.tt;eps = Stat.vpro3.tke;
eps(eps == 0) = -inf;
eps(isnan(eps)) = -inf;
imagesc(tt,z,eps)
caxis([1E-5 3E-4])
datetickzoom('x','HH:MM','keepticks','keeplimits')
set(gca,'ydir','normal')


cb2 = colorbar;
freezeColors
cdh = get(cb2,'children');
map = get(cdh,'cdata');
set(cdh,'cdata',map(2:end));
cb2 = cbfreeze(cb2);
set(cb2,'YTick',1E-4:1E-4:3E-4)


%%%epsilon subsections
tt = Stat.vpro1.tt;
stop2 = datenum(2015,03,06,14,12,00);
ind = find(tt >= start & tt <= stop2);
z = linspace(5,20,15);z = fliplr(z);

sp(7) = subplot(337);
tt = Stat.vpro1.tt(ind);eps = Stat.vpro1.tke(1:15,ind);
eps(eps == 0) = -inf;
eps(isnan(eps)) = -inf;
imagesc(tt,z,eps)
caxis([1E-5 3E-4])
freezeColors
datetickzoom('x','HH:MM','keepticks','keeplimits')
set(gca,'ydir','normal')

sp(8) = subplot(338);
tt = Stat.vpro2.tt(ind);eps = Stat.vpro2.tke(1:15,ind);
eps(eps == 0) = -inf;
eps(isnan(eps)) = -inf;
imagesc(tt,z,eps)
caxis([1E-5 3E-4])
freezeColors
datetickzoom('x','HH:MM','keepticks','keeplimits')
set(gca,'ydir','normal')

sp(9) = subplot(339);
tt = Stat.vpro3.tt(ind);eps = Stat.vpro3.tke(1:15,ind);
eps(eps == 0) = -inf;
eps(isnan(eps)) = -inf;
imagesc(tt,z,eps)
caxis([1E-5 3E-4])
datetickzoom('x','HH:MM','keepticks','keeplimits')
set(gca,'ydir','normal')

cb3 = colorbar;
freezeColors
cdh = get(cb3,'children');
map = get(cdh,'cdata');
set(cdh,'cdata',map(2:end));
cb3 = cbfreeze(cb3);
set(cb3,'YTick',1E-4:1E-4:3E-4)


%%%set positioning, axes adjustments
set([sp(1) sp(2) sp(3)...
    sp(4) sp(5) sp(6)],...
    'XTick',start:2*twm:stop,...
    'Xlim',[start stop],...
    'Ylim',[-10 20],...
    'Box','on',...
    'LineWidth',1.5,...
    'FontSize',14,...
    'FontName','Arial',...
    'TickDir','out')
set([sp(7) sp(8) sp(9)],...
    'XTick',start:ten:stop2,...
    'Xlim',[start stop2],...
    'Ylim',[5 15],...
    'Box','on',...
    'LineWidth',1.5,...
    'FontSize',14,...
    'FontName','Arial',...
    'TickDir','out')

set([sp(1) sp(2) sp(3)],...
    'XTickLabel',[])
set([sp(2) sp(3) sp(5) sp(6)],...
    'YTickLabel',[])
set(sp(1),'position',[0.09 0.72 0.22 0.24])
set(sp(2),'position',[0.36 0.72 0.22 0.24])
set(sp(3),'position',[0.63 0.72 0.22 0.24])
set(sp(4),'position',[0.09 0.43 0.22 0.24])
set(sp(5),'position',[0.36 0.43 0.22 0.24])
set(sp(6),'position',[0.63 0.43 0.22 0.24])
set(sp(7),'position',[0.09 0.1 0.22 0.24])
set(sp(8),'position',[0.36 0.1 0.22 0.24])
set(sp(9),'position',[0.63 0.1 0.22 0.24])

set(cb1,'position',[0.88 0.72 0.025 0.24])
set(cb2,'position',[0.88 0.43 0.025 0.24])
set(cb3,'position',[0.88 0.1 0.025 0.24])

set([cb1 cb2 cb3],...
    'LineWidth',1.5,...
    'TickDir','out',...
    'FontSize',14,...
    'FontName','Arial')
ylabel(cb1,'Onshore Velocity (ms^-^1)','FontSize',18)
ylabel(cb2,'\epsilon (Wkg^-^1)','FontSize',18)
ylabel(cb3,'\epsilon (Wkg^-^1)','FontSize',18)
xlabel(sp(5),'Time on March 6th, 2015','FontSize',18)
xlabel(sp(8),'Time on March 6th, 2015','FontSize',18)
ylabel(sp(1),'HAB (mm)','FontSize',18)
ylabel(sp(4),'HAB (mm)','FontSize',18)
ylabel(sp(7),'HAB (mm)','FontSize',18)
fdir = 'd:\Projects\Mekong_W2015\Figures\Paper1\CurrentDirections\';
export_fig([fdir 'AveVelsTKE_TS_bw'],'-pdf','-nocrop')