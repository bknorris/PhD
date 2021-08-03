clear
load('C:\Users\bkn5\Projects\Mekong_W2015\Data\Aquadopp\FSS\AD5116_12March2015.mat')

%time average vertical velocity:
start = datenum(2015,03,10,15,00,00);stop = datenum(2015,03,10,16,40,00);step = datenum(0,0,0,0,30,0);
heading = 110; %cross shore angle
idx = find(aqdp.datenum >= start & aqdp.datenum <= stop);
h = 51.3; %inst hab (cm)
ph = 41.4; %pneumatophore height (cm)

b1 = aqdp.u(idx,1:7);
b2 = aqdp.v(idx,1:7);
b3 = aqdp.w(idx,1:7);
t = aqdp.datenum(idx,:);
p = aqdp.pressure(idx,:);

%rotate to the cross-shore direction
rot = heading*pi/180;
u = b1.*(ones(size(b1))*cos(rot)) + ...
         b2.*(ones(size(b2))*sin(rot));
v = -b1.*(ones(size(b1))*sin(rot)) + ...
        b2.*(ones(size(b2))*cos(rot));
w = b3;
rb = aqdp.rangebins;rbcm = h-rb(1:7)*100;
sr = 8;
intv = 2; %minutes
avt = intv*sr*60;
ind = [1 avt:avt:length(w)];

f1 = figure(1);
set(f1,'PaperOrientation','portrait',...
    'position',[200 200   1000   600]);
p(1) = subplot(131);
for i = 1:length(ind)-1
    Uz = u(ind(i):ind(i+1),:);
    U(i,:) = nanmean(Uz);
    T(i,:) = t(ind(i));
end
% c = jet(length(ind));
imagesc(T,rbcm,U'*100), hold on
phx = linspace(T(1),T(end),length(T));
plot(phx,ones(length(phx))*ph,'--k','LineWidth',1.5)
set(gca,'Ydir','normal','YGrid','on','YLim',[28 46],'xlim',[T(1) T(end)],'xtick',T(1):step:T(end))
caxis([-5 5])
datetickzoom('x','HH:MM','keepticks','keeplimits')
xlabel('Velocity (cm/s)')
ylabel('Height Above Bed (cm)')
title('Along-Shore')

p(2) = subplot(132);
for i = 1:length(ind)-1
    Vz = v(ind(i):ind(i+1),:);
    V(i,:) = nanmean(Vz);
end
imagesc(T,rbcm,V'*100), hold on
phx = linspace(T(1),T(end),length(T));
plot(phx,ones(length(phx))*ph,'--k','LineWidth',1.5)
set(gca,'Ydir','normal','YGrid','on','YTickLabel',[],'YLim',[28 46],'xlim',[T(1) T(end)],'xtick',T(1):step:T(end))
caxis([-5 5])
datetickzoom('x','HH:MM','keepticks','keeplimits')
xlabel('Velocity (cm/s)')
% ylabel('Height Above Bed (cm)')
title('Cross-Shore')

p(3) = subplot(133);
for i = 1:length(ind)-1
    Wz = w(ind(i):ind(i+1),:);
    W(i,:) = nanmean(Wz);
end
imagesc(T,rbcm,W'*100), hold on
phx = linspace(T(1),T(end),length(T));
plot(phx,ones(length(phx))*ph,'--k','LineWidth',1.5)
set(gca,'Ydir','normal','YGrid','on','YTickLabel',[],'YLim',[28 46],'xlim',[T(1) T(end)],'xtick',T(1):step:T(end))
caxis([-5 5])
datetickzoom('x','HH:MM','keepticks','keeplimits')
xlabel('Velocity (cm/s)')
% ylabel('Height Above Bed (cm)')
title('Vertical')

set(p(1),'position',[0.1 0.22 0.25 0.7])
set(p(2),'position',[0.38 0.22 0.25 0.7])
set(p(3),'position',[0.66 0.22 0.25 0.7])
cb = colorbar('south');
set(cb,'position',[0.25 0.02, 0.54 0.05]),xlabel(cb,'cm/s')
suptitle(['AD5116 - Above Bare-Bed: ' datestr(start,'dd/mm HH:MM')])

f2 = figure(2);
subplot(131)
Ubar = mean(U,1);
stdv = std(U,0,1);
rbcm = h-rb(1:7)*100;
p(1) = plot(Ubar*100,rbcm,...
    'Color','k','Marker','^','MarkerFaceColor','k','LineWidth',1.5);hold on
p(2) = plot((Ubar-stdv)*100,rbcm,...
    'Color',[0.6 0.6 0.6],'LineStyle','--','LineWidth',1.5);
p(3) = plot((Ubar+stdv)*100,rbcm,...
    'Color',[0.6 0.6 0.6],'LineStyle','--','LineWidth',1.5);
set(gca,'Xlim',[-20 10])
title('Along-Shore')
xlabel('Velocity (cm/s)')
ylabel('Height Above Bed (cm)')

subplot(132)
Vbar = mean(V,1);
stdv = std(V,0,1);
p(1) = plot(Vbar*100,rbcm,...
    'Color','k','Marker','o','MarkerFaceColor','k','LineWidth',1.5);hold on
p(2) = plot((Vbar-stdv)*100,rbcm,...
    'Color',[0.6 0.6 0.6],'LineStyle','--','LineWidth',1.5);
p(3) = plot((Vbar+stdv)*100,rbcm,...
    'Color',[0.6 0.6 0.6],'LineStyle','--','LineWidth',1.5);
set(gca,'Xlim',[-10 10])
title('Cross-Shore')
xlabel('Velocity (cm/s)')
ylabel('Height Above Bed (cm)')

subplot(133)
Wbar = mean(W,1);
stdv = std(W,0,1);
p(1) = plot(Wbar*100,rbcm,...
    'Color','k','Marker','d','MarkerFaceColor','k','LineWidth',1.5);hold on
p(2) = plot((Wbar-stdv)*100,rbcm,...
    'Color',[0.6 0.6 0.6],'LineStyle','--','LineWidth',1.5);
p(3) = plot((Wbar+stdv)*100,rbcm,...
    'Color',[0.6 0.6 0.6],'LineStyle','--','LineWidth',1.5);
set(gca,'Xlim',[-5 5])
title('Vertical')
xlabel('Velocity (cm/s)')
ylabel('Height Above Bed (cm)')
suptitle('AD5116 - Above Bare-Bed')

clear U V W T
%%%%%%%

load('C:\Users\bkn5\Projects\Mekong_W2015\Data\Aquadopp\FSS\AD5117_12March2015.mat')
idx = find(aqdp.datenum >= start & aqdp.datenum <= stop);
h = 54.8; %inst hab (cm)
b1 = aqdp.u(idx,1:7);
b2 = aqdp.v(idx,1:7);
b3 = aqdp.w(idx,1:7);
t = aqdp.datenum(idx,:);
p = aqdp.pressure(idx,:);

%rotate to the cross-shore direction
rot = heading*pi/180;
u = b1.*(ones(size(b1))*cos(rot)) + ...
         b2.*(ones(size(b2))*sin(rot));
v = -b1.*(ones(size(b1))*sin(rot)) + ...
        b2.*(ones(size(b2))*cos(rot));
w = b3;
rb = aqdp.rangebins;rbcm = h-rb(1:7)*100;
sr = 8;
avt = intv*sr*60;
ind = [1 avt:avt:length(w)];

f3 = figure(3);
set(f3,'PaperOrientation','portrait',...
    'position',[200 200   1000   600]);
p(1) = subplot(131);
c = jet(length(ind));
for i = 1:length(ind)-1
    Uz = u(ind(i):ind(i+1),:);
    U(i,:) = nanmean(Uz);
    T(i,:) = t(ind(i));
end
imagesc(T,rbcm,U'*100), hold on
phx = linspace(T(1),T(end),length(T));
plot(phx,ones(length(phx))*ph,'--k','LineWidth',1.5)
set(gca,'Ydir','normal','YGrid','on','YLim',[32 50],'xlim',[T(1) T(end)],'xtick',T(1):step:T(end))
caxis([-5 5])
datetickzoom('x','HH:MM','keepticks','keeplimits')
xlabel('Velocity (cm/s)')
ylabel('Height Above Bed (cm)')
title('Along-Shore')

p(2) = subplot(132);
for i = 1:length(ind)-1
    Vz = v(ind(i):ind(i+1),:);
    V(i,:) = nanmean(Vz);
end
imagesc(T,rbcm,V'*100), hold on
phx = linspace(T(1),T(end),length(T));
plot(phx,ones(length(phx))*ph,'--k','LineWidth',1.5)
set(gca,'Ydir','normal','YGrid','on','YTickLabel',[],'YLim',[32 50],'xlim',[T(1) T(end)],'xtick',T(1):step:T(end))
caxis([-5 5])
datetickzoom('x','HH:MM','keepticks','keeplimits')
xlabel('Velocity (cm/s)')
% ylabel('Height Above Bed (cm)')
title('Cross-Shore')

p(3) = subplot(133);
for i = 1:length(ind)-1
    Wz = w(ind(i):ind(i+1),:);
    W(i,:) = nanmean(Wz);
end
imagesc(T,rbcm,W'*100), hold on
phx = linspace(T(1),T(end),length(T));
plot(phx,ones(length(phx))*ph,'--k','LineWidth',1.5)
set(gca,'Ydir','normal','YGrid','on','YTickLabel',[],'YLim',[32 50],'xlim',[T(1) T(end)],'xtick',T(1):step:T(end))
caxis([-5 5])
datetickzoom('x','HH:MM','keepticks','keeplimits')
xlabel('Velocity (cm/s)')
% ylabel('Height Above Bed (cm)')
title('Vertical')

set(p(1),'position',[0.1 0.22 0.25 0.7])
set(p(2),'position',[0.38 0.22 0.25 0.7])
set(p(3),'position',[0.66 0.22 0.25 0.7])
cb = colorbar('south');
set(cb,'position',[0.25 0.02, 0.54 0.05]),xlabel(cb,'cm/s')
suptitle(['AD5117 - Above Pneumatophores: ' datestr(start,'dd/mm HH:MM')])

f4 = figure(4);
subplot(131)
Ubar = mean(U,1);
stdv = std(U,0,1);
rbcm = h-rb(1:7)*100;
p(1) = plot(Ubar*100,rbcm,...
    'Color','k','Marker','^','MarkerFaceColor','k','LineWidth',1.5);hold on
p(2) = plot((Ubar-stdv)*100,rbcm,...
    'Color',[0.6 0.6 0.6],'LineStyle','--','LineWidth',1.5);
p(3) = plot((Ubar+stdv)*100,rbcm,...
    'Color',[0.6 0.6 0.6],'LineStyle','--','LineWidth',1.5);
set(gca,'Xlim',[-20 10])
title('Along-Shore')
xlabel('Velocity (cm/s)')
ylabel('Height Above Bed (cm)')

subplot(132)
Vbar = mean(V,1);
stdv = std(V,0,1);
p(1) = plot(Vbar*100,rbcm,...
    'Color','k','Marker','o','MarkerFaceColor','k','LineWidth',1.5);hold on
p(2) = plot((Vbar-stdv)*100,rbcm,...
    'Color',[0.6 0.6 0.6],'LineStyle','--','LineWidth',1.5);
p(3) = plot((Vbar+stdv)*100,rbcm,...
    'Color',[0.6 0.6 0.6],'LineStyle','--','LineWidth',1.5);
set(gca,'Xlim',[-10 10])
title('Cross-Shore')
xlabel('Velocity (cm/s)')
ylabel('Height Above Bed (cm)')

subplot(133)
Wbar = mean(W,1);
stdv = std(W,0,1);
p(1) = plot(Wbar*100,rbcm,...
    'Color','k','Marker','d','MarkerFaceColor','k','LineWidth',1.5);hold on
p(2) = plot((Wbar-stdv)*100,rbcm,...
    'Color',[0.6 0.6 0.6],'LineStyle','--','LineWidth',1.5);
p(3) = plot((Wbar+stdv)*100,rbcm,...
    'Color',[0.6 0.6 0.6],'LineStyle','--','LineWidth',1.5);
set(gca,'Xlim',[-5 5])
title('Vertical')
xlabel('Velocity (cm/s)')
ylabel('Height Above Bed (cm)')
suptitle('AD5117 - Above Pneumatophores')

% clear U V W T
%%%%%%%

% load('C:\Users\bkn5\Projects\Mekong_W2015\Data\Aquadopp\FSS\ADLR4_9March2015_f.mat')
% idx = find(aqdp.datenum >= start & aqdp.datenum <= stop);
% h = 61; %inst hab (cm)
% b1 = aqdp.u(idx,:);
% b2 = aqdp.v(idx,:);
% b3 = aqdp.w(idx,:);
% t = aqdp.datenum(idx,:);
% p = aqdp.pressure(idx,:);
% 
% %rotate to the cross-shore direction
% rot = heading*pi/180;
% u = b1.*(ones(size(b1))*cos(rot)) + ...
%          b2.*(ones(size(b2))*sin(rot));
% v = -b1.*(ones(size(b1))*sin(rot)) + ...
%         b2.*(ones(size(b2))*cos(rot));
% w = b3;
% rb = aqdp.rangebins;
% sr = 1;
% avt = intv*sr*60;
% ind = [1 avt:avt:length(w)];
% 
% f5 = figure(5);
% set(f5,'PaperOrientation','portrait',...
%     'position',[200 200   1000   600]);
% p(1) = subplot(131);
% c = jet(length(ind));
% for i = 1:length(ind)-1
%     Uz = u(ind(i):ind(i+1),:);
%     U(i,:) = nanmean(Uz);
%     T(i,:) = t(ind(i));
% end
% imagesc(T,rbcm,U'*100)
% set(gca,'Ydir','normal','YGrid','on')
% caxis([-0.01 0.01])
% datetickzoom('x','HH:MM','keepticks','keeplimits')
% xlabel('Velocity (cm/s)')
% ylabel('Height Above Bed (cm)')
% title('Along-Shore')
% 
% p(2) = subplot(132);
% for i = 1:length(ind)-1
%     Vz = v(ind(i):ind(i+1),:);
%     V(i,:) = nanmean(Vz);
% end
% imagesc(T,rbcm,V'*100)
% set(gca,'Ydir','normal','YGrid','on','YTickLabel',[])
% caxis([-0.01 0.01])
% datetickzoom('x','HH:MM','keepticks','keeplimits')
% xlabel('Velocity (cm/s)')
% % ylabel('Height Above Bed (cm)')
% title('Cross-Shore')
% 
% p(3) = subplot(133);
% for i = 1:length(ind)-1
%     Wz = w(ind(i):ind(i+1),:);
%     W(i,:) = nanmean(Wz);
% end
% imagesc(T,rbcm,W'*100)
% set(gca,'Ydir','normal','YGrid','on','YTickLabel',[])
% caxis([-0.01 0.01])
% datetickzoom('x','HH:MM','keepticks','keeplimits')
% xlabel('Velocity (cm/s)')
% % ylabel('Height Above Bed (cm)')
% title('Vertical')
% 
% set(p(1),'position',[0.1 0.22 0.25 0.7])
% set(p(2),'position',[0.38 0.22 0.25 0.7])
% set(p(3),'position',[0.66 0.22 0.25 0.7])
% cb = colorbar('south');
% set(cb,'position',[0.25 0.02, 0.54 0.05]),xlabel(cb,'cm/s')
% 
% f6 = figure(6);
% subplot(131)
% Ubar = mean(U,1);
% stdv = std(U,0,1);
% rbcm = rbcm;
% p(1) = plot(Ubar*100,rbcm,...
%     'Color','k','Marker','^','MarkerFaceColor','k','LineWidth',1.5);hold on
% p(2) = plot((Ubar-stdv)*100,rbcm,...
%     'Color',[0.6 0.6 0.6],'LineStyle','--','LineWidth',1.5);
% p(3) = plot((Ubar+stdv)*100,rbcm,...
%     'Color',[0.6 0.6 0.6],'LineStyle','--','LineWidth',1.5);
% set(gca,'Xlim',[-0.02 0])
% title('Along-Shore')
% xlabel('Velocity (cm/s)')
% ylabel('Height Above Bed (cm)')
% 
% subplot(132)
% Vbar = mean(V,1);
% stdv = std(V,0,1);
% p(1) = plot(Vbar*100,rbcm,...
%     'Color','k','Marker','o','MarkerFaceColor','k','LineWidth',1.5);hold on
% p(2) = plot((Vbar-stdv)*100,rbcm,...
%     'Color',[0.6 0.6 0.6],'LineStyle','--','LineWidth',1.5);
% p(3) = plot((Vbar+stdv)*100,rbcm,...
%     'Color',[0.6 0.6 0.6],'LineStyle','--','LineWidth',1.5);
% set(gca,'Xlim',[-3 1])
% title('Cross-Shore')
% xlabel('Velocity (cm/s)')
% ylabel('Height Above Bed (cm)')
% 
% subplot(133)
% Wbar = mean(W,1);
% stdv = std(W,0,1);
% p(1) = plot(Wbar*100,rbcm,...
%     'Color','k','Marker','d','MarkerFaceColor','k','LineWidth',1.5);hold on
% p(2) = plot((Wbar-stdv)*100,rbcm,...
%     'Color',[0.6 0.6 0.6],'LineStyle','--','LineWidth',1.5);
% p(3) = plot((Wbar+stdv)*100,rbcm,...
%     'Color',[0.6 0.6 0.6],'LineStyle','--','LineWidth',1.5);
% set(gca,'Xlim',[-1 3])
% title('Vertical')
% xlabel('Velocity (cm/s)')
% ylabel('Height Above Bed (cm)')