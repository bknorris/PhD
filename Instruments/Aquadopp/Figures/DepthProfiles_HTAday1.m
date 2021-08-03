clear
load('C:\Users\bkn5\Projects\Mekong_W2015\Data\Aquadopp\FSS\AD5116_9March2015.mat')

%time average vertical velocity:
start = datenum(2015,03,07,14,15,00);stop = datenum(2015,03,07,17,10,00);step = datenum(0,0,0,0,30,0);
heading = 110; %cross shore angle
idx = find(aqdp.datenum >= start & aqdp.datenum <= stop);
h = 63.2; %inst hab (cm)
ph = 43.5; %pneumatophore height (cm)

b1 = aqdp.u(idx,5:end);
b2 = aqdp.v(idx,5:end);
b3 = aqdp.w(idx,5:end);
t = aqdp.datenum(idx,:);
p = aqdp.pressure(idx,:);

%rotate to the cross-shore direction
rot = heading*pi/180;
u = b1.*(ones(size(b1))*cos(rot)) + ...
         b2.*(ones(size(b2))*sin(rot));
v = -b1.*(ones(size(b1))*sin(rot)) + ...
        b2.*(ones(size(b2))*cos(rot));
w = b3;
rb = aqdp.rangebins;rbcm = h-rb(5:end)*100;
dat.a6.rbcm = rbcm;
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
dat.a6.U = U;
% c = jet(length(ind));
imagesc(T,rbcm,U'*100)
set(gca,'Ydir','normal','YGrid','on','YLim',[16 46],'xlim',[T(1) T(end)],'xtick',T(1):step:T(end))
caxis([-5 5])
datetickzoom('x','HH:MM','keepticks','keeplimits')
xlabel('Time on 07/03/16')
ylabel('Height Above Bed (cm)')
title('Along-Shore')

p(2) = subplot(132);
for i = 1:length(ind)-1
    Vz = v(ind(i):ind(i+1),:);
    V(i,:) = nanmean(Vz);
end
dat.a6.V = V;
imagesc(T,rbcm,V'*100)
set(gca,'Ydir','normal','YGrid','on','YTickLabel',[],'YLim',[16 46],'xlim',[T(1) T(end)],'xtick',T(1):step:T(end))
caxis([-5 5])
datetickzoom('x','HH:MM','keepticks','keeplimits')
xlabel('Time on 07/03/16')
% ylabel('Height Above Bed (cm)')
title('Cross-Shore')

p(3) = subplot(133);
for i = 1:length(ind)-1
    Wz = w(ind(i):ind(i+1),:);
    W(i,:) = nanmean(Wz);
end
dat.a6.W = W;
imagesc(T,rbcm,W'*100)
set(gca,'Ydir','normal','YGrid','on','YTickLabel',[],'YLim',[16 46],'xlim',[T(1) T(end)],'xtick',T(1):step:T(end))
caxis([-5 5])
datetickzoom('x','HH:MM','keepticks','keeplimits')
xlabel('Time on 07/03/16')
% ylabel('Height Above Bed (cm)')
title('Vertical')

set(p(1),'position',[0.1 0.22 0.25 0.7])
set(p(2),'position',[0.38 0.22 0.25 0.7])
set(p(3),'position',[0.66 0.22 0.25 0.7])
cb = colorbar('south');
set(cb,'position',[0.25 0.02, 0.54 0.05]),xlabel(cb,'cm/s')
suptitle(['AD5116 - Above Bare-Bed: ' datestr(start,'dd/mm HH:MM')])

f2 = figure(2);
phx = linspace(-50,50,10);
subplot(131)
Ubar = mean(U,1);
stdv = std(U,0,1);
rbcm = h-rb(5:end)*100;
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
set(gca,'Xlim',[-2 2])
title('Vertical')
xlabel('Velocity (cm/s)')
ylabel('Height Above Bed (cm)')
suptitle('AD5116 - Above Bare-Bed')

%Calculate RMS Velocities
f3 = figure(3);
set(f3,'PaperOrientation','portrait',...
    'position',[200 200   1000   600]);
p(1) = subplot(131);
for i = 1:length(ind)-1
    Uz = u(ind(i):ind(i+1),:);
    Uz(isnan(Uz)) = 0; %remove NaNs- this will not affect the RMS calc
    Urms(i,:) = sqrt((sum(Uz.^2)/numel(Uz)));
    T(i,:) = t(ind(i));
end
dat.a6.Urms = Urms;
% c = jet(length(ind));
imagesc(T,rbcm,Urms'*100)
set(gca,'Ydir','normal','YGrid','on','YLim',[16 46],'xlim',[T(1) T(end)],'xtick',T(1):step:T(end))
caxis([0 6])
datetickzoom('x','HH:MM','keepticks','keeplimits')
xlabel('Time on 07/03/16')
ylabel('Height Above Bed (cm)')
title('Along-Shore RMS Velocities')

p(2) = subplot(132);
for i = 1:length(ind)-1
    Vz = v(ind(i):ind(i+1),:);
    Vz(isnan(Vz)) = 0; %remove NaNs- this will not affect the RMS calc
    Vrms(i,:) = sqrt((sum(Vz.^2)/numel(Vz)));
    T(i,:) = t(ind(i));
end
dat.a6.Vrms = Vrms;
% c = jet(length(ind));
imagesc(T,rbcm,Vrms'*100)
set(gca,'Ydir','normal','YGrid','on','YTickLabel',[],'YLim',[16 46],'xlim',[T(1) T(end)],'xtick',T(1):step:T(end))
caxis([0 6])
datetickzoom('x','HH:MM','keepticks','keeplimits')
xlabel('Time on 07/03/16')
title('Cross-Shore RMS Velocities')

p(3) = subplot(133);
for i = 1:length(ind)-1
    Wz = w(ind(i):ind(i+1),:);
    Wz(isnan(Wz)) = 0; %remove NaNs- this will not affect the RMS calc
    Wrms(i,:) = sqrt((sum(Wz.^2)/numel(Wz)));
    T(i,:) = t(ind(i));
end
% c = jet(length(ind));
dat.a6.Wrms = Wrms;
imagesc(T,rbcm,Wrms'*100)
set(gca,'Ydir','normal','YGrid','on','YTickLabel',[],'YLim',[16 46],'xlim',[T(1) T(end)],'xtick',T(1):step:T(end))
caxis([0 6])
datetickzoom('x','HH:MM','keepticks','keeplimits')
xlabel('Time on 07/03/16')
title('Vertical RMS Velocities')
set(p(1),'position',[0.1 0.22 0.25 0.7])
set(p(2),'position',[0.38 0.22 0.25 0.7])
set(p(3),'position',[0.66 0.22 0.25 0.7])
cb = colorbar('south');
set(cb,'position',[0.25 0.02, 0.54 0.05]),xlabel(cb,'cm/s')
suptitle(['AD5116 - Above Bare-Bed: ' datestr(start,'dd/mm HH:MM')])

f4 = figure(4);
phx = linspace(-50,50,10);
subplot(131)
Unorm = U./Urms;
Ubar = mean(Unorm,1); 
stdv = std(Unorm,0,1);
rbcm = h-rb(5:end)*100;
p(1) = plot(Ubar,rbcm,...
    'Color','k','Marker','^','MarkerFaceColor','k','LineWidth',1.5);hold on
p(2) = plot((Ubar-stdv),rbcm,...
    'Color',[0.6 0.6 0.6],'LineStyle','--','LineWidth',1.5);
p(3) = plot((Ubar+stdv),rbcm,...
    'Color',[0.6 0.6 0.6],'LineStyle','--','LineWidth',1.5);
set(gca,'Xlim',[-1 2])
title('Along-Shore')
xlabel('U/U_r_m_s')
ylabel('Height Above Bed (cm)')

subplot(132)
Vnorm = V./Vrms;
Vbar = mean(Vnorm,1);
stdv = std(Vnorm,0,1);
p(1) = plot(Vbar,rbcm,...
    'Color','k','Marker','o','MarkerFaceColor','k','LineWidth',1.5);hold on
p(2) = plot((Vbar-stdv),rbcm,...
    'Color',[0.6 0.6 0.6],'LineStyle','--','LineWidth',1.5);
p(3) = plot((Vbar+stdv),rbcm,...
    'Color',[0.6 0.6 0.6],'LineStyle','--','LineWidth',1.5);
set(gca,'Xlim',[0 3])
title('Cross-Shore')
xlabel('V/V_r_m_s')
ylabel('Height Above Bed (cm)')

subplot(133)
Wnorm = W./Wrms;
Wbar = mean(Wnorm,1);
stdv = std(Wnorm,0,1);
p(1) = plot(Wbar,rbcm,...
    'Color','k','Marker','d','MarkerFaceColor','k','LineWidth',1.5);hold on
p(2) = plot((Wbar-stdv),rbcm,...
    'Color',[0.6 0.6 0.6],'LineStyle','--','LineWidth',1.5);
p(3) = plot((Wbar+stdv),rbcm,...
    'Color',[0.6 0.6 0.6],'LineStyle','--','LineWidth',1.5);
set(gca,'Xlim',[-1 0])
title('Vertical')
xlabel('W/W_r_m_s')
ylabel('Height Above Bed (cm)')
suptitle('AD5116 - Above Bare-Bed')

clear U Urms V Vrms W Wrms T
%%%%%%%

load('C:\Users\bkn5\Projects\Mekong_W2015\Data\Aquadopp\FSS\AD5117_9March2015.mat')
idx = find(aqdp.datenum >= start & aqdp.datenum <= stop);
h = 68.7; %inst hab (cm)
b1 = aqdp.u(idx,7:end);
b2 = aqdp.v(idx,7:end);
b3 = aqdp.w(idx,7:end);
t = aqdp.datenum(idx,:);
p = aqdp.pressure(idx,:);

%rotate to the cross-shore direction
rot = heading*pi/180;
u = b1.*(ones(size(b1))*cos(rot)) + ...
         b2.*(ones(size(b2))*sin(rot));
v = -b1.*(ones(size(b1))*sin(rot)) + ...
        b2.*(ones(size(b2))*cos(rot));
w = b3;
rb = aqdp.rangebins;rbcm = h-rb(7:end)*100;
dat.a7.rbcm = rbcm;
sr = 8;
avt = intv*sr*60;
ind = [1 avt:avt:length(w)];

f5 = figure(5);
set(f5,'PaperOrientation','portrait',...
    'position',[200 200   1000   600]);
p(1) = subplot(131);
c = jet(length(ind));
for i = 1:length(ind)-1
    Uz = u(ind(i):ind(i+1),:);
    U(i,:) = nanmean(Uz);
    T(i,:) = t(ind(i));
end
dat.a7.U = U;
imagesc(T,rbcm,U'*100), hold on
phx = linspace(T(1),T(end),length(T));
plot(phx,ones(length(phx))*ph,'--k','LineWidth',1.5)
set(gca,'Ydir','normal','YGrid','on','YLim',[22 45],'xlim',[T(1) T(end)],'xtick',T(1):step:T(end))
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
dat.a7.V = V;
imagesc(T,rbcm,V'*100), hold on
phx = linspace(T(1),T(end),length(T));
plot(phx,ones(length(phx))*ph,'--k','LineWidth',1.5)
set(gca,'Ydir','normal','YGrid','on','YTickLabel',[],'YLim',[22 45],'xlim',[T(1) T(end)],'xtick',T(1):step:T(end))
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
dat.a7.W = W;
imagesc(T,rbcm,W'*100), hold on
phx = linspace(T(1),T(end),length(T));
plot(phx,ones(length(phx))*ph,'--k','LineWidth',1.5)
set(gca,'Ydir','normal','YGrid','on','YTickLabel',[],'YLim',[22 45],'xlim',[T(1) T(end)],'xtick',T(1):step:T(end))
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

f6 = figure(6);
phx = linspace(-50,50,10);
subplot(131)
Ubar = mean(U,1);
stdv = std(U,0,1);
rbcm = h-rb(7:end)*100;
p(1) = plot(Ubar*100,rbcm,...
    'Color','k','Marker','^','MarkerFaceColor','k','LineWidth',1.5);hold on
p(2) = plot((Ubar-stdv)*100,rbcm,...
    'Color',[0.6 0.6 0.6],'LineStyle','--','LineWidth',1.5);
p(3) = plot((Ubar+stdv)*100,rbcm,...
    'Color',[0.6 0.6 0.6],'LineStyle','--','LineWidth',1.5);
plot(phx,ones(length(phx))*ph,'--k','LineWidth',1.5)
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
plot(phx,ones(length(phx))*ph,'--k','LineWidth',1.5)
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
plot(phx,ones(length(phx))*ph,'--k','LineWidth',1.5)
set(gca,'Xlim',[-2 2])
title('Vertical')
xlabel('Velocity (cm/s)')
ylabel('Height Above Bed (cm)')
suptitle('AD5117 - Above Pneumatophores')

%Calculate RMS Velocities
f7 = figure(7);
set(f7,'PaperOrientation','portrait',...
    'position',[200 200   1000   600]);
p(1) = subplot(131);
for i = 1:length(ind)-1
    Uz = u(ind(i):ind(i+1),:);
    Uz(isnan(Uz)) = 0; %remove NaNs- this will not affect the RMS calc
    Urms(i,:) = sqrt((sum(Uz.^2)/numel(Uz)));
    T(i,:) = t(ind(i));
end
dat.a7.Urms = Urms;
% c = jet(length(ind));
imagesc(T,rbcm,Urms'*100), hold on
phx = linspace(T(1),T(end),length(T));
plot(phx,ones(length(phx))*ph,'--k','LineWidth',1.5)
set(gca,'Ydir','normal','YGrid','on','YLim',[22 45],'xlim',[T(1) T(end)],'xtick',T(1):step:T(end))
caxis([0 6])
datetickzoom('x','HH:MM','keepticks','keeplimits')
xlabel('Time on 07/03/16')
ylabel('Height Above Bed (cm)')
title('Along-Shore RMS Velocities')

p(2) = subplot(132);
for i = 1:length(ind)-1
    Vz = v(ind(i):ind(i+1),:);
    Vz(isnan(Vz)) = 0; %remove NaNs- this will not affect the RMS calc
    Vrms(i,:) = sqrt((sum(Vz.^2)/numel(Vz)));
    T(i,:) = t(ind(i));
end
dat.a7.Vrms = Vrms;
% c = jet(length(ind));
imagesc(T,rbcm,Vrms'*100), hold on
phx = linspace(T(1),T(end),length(T));
plot(phx,ones(length(phx))*ph,'--k','LineWidth',1.5)
set(gca,'Ydir','normal','YGrid','on','YTickLabel',[],'YLim',[22 45],'xlim',[T(1) T(end)],'xtick',T(1):step:T(end))
caxis([0 6])
datetickzoom('x','HH:MM','keepticks','keeplimits')
xlabel('Time on 07/03/16')
title('Cross-Shore RMS Velocities')

p(3) = subplot(133);
for i = 1:length(ind)-1
    Wz = w(ind(i):ind(i+1),:);
    Wz(isnan(Wz)) = 0; %remove NaNs- this will not affect the RMS calc
    Wrms(i,:) = sqrt((sum(Wz.^2)/numel(Wz)));
    T(i,:) = t(ind(i));
end
dat.a7.Wrms = Wrms;
% c = jet(length(ind));
imagesc(T,rbcm,Wrms'*100), hold on
phx = linspace(T(1),T(end),length(T));
plot(phx,ones(length(phx))*ph,'--k','LineWidth',1.5)
set(gca,'Ydir','normal','YGrid','on','YTickLabel',[],'YLim',[22 45],'xlim',[T(1) T(end)],'xtick',T(1):step:T(end))
caxis([0 6])
datetickzoom('x','HH:MM','keepticks','keeplimits')
xlabel('Time on 07/03/16')
title('Vertical RMS Velocities')
set(p(1),'position',[0.1 0.22 0.25 0.7])
set(p(2),'position',[0.38 0.22 0.25 0.7])
set(p(3),'position',[0.66 0.22 0.25 0.7])
cb = colorbar('south');
set(cb,'position',[0.25 0.02, 0.54 0.05]),xlabel(cb,'cm/s')
suptitle(['AD5117 - Above Pneumatophores: ' datestr(start,'dd/mm HH:MM')])

f8 = figure(8);
phx = linspace(-50,50,10);
subplot(131)
Unorm = U./Urms;
Ubar = mean(Unorm,1); 
stdv = std(Unorm,0,1);
p(1) = plot(Ubar,rbcm,...
    'Color','k','Marker','^','MarkerFaceColor','k','LineWidth',1.5);hold on
p(2) = plot((Ubar-stdv),rbcm,...
    'Color',[0.6 0.6 0.6],'LineStyle','--','LineWidth',1.5);
p(3) = plot((Ubar+stdv),rbcm,...
    'Color',[0.6 0.6 0.6],'LineStyle','--','LineWidth',1.5);
plot(phx,ones(length(phx))*ph,'--k','LineWidth',1.5)
set(gca,'Xlim',[-1 2])
title('Along-Shore')
xlabel('U/U_r_m_s')
ylabel('Height Above Bed (cm)')

subplot(132)
Vnorm = V./Vrms;
Vbar = mean(Vnorm,1);
stdv = std(Vnorm,0,1);
p(1) = plot(Vbar,rbcm,...
    'Color','k','Marker','o','MarkerFaceColor','k','LineWidth',1.5);hold on
p(2) = plot((Vbar-stdv),rbcm,...
    'Color',[0.6 0.6 0.6],'LineStyle','--','LineWidth',1.5);
p(3) = plot((Vbar+stdv),rbcm,...
    'Color',[0.6 0.6 0.6],'LineStyle','--','LineWidth',1.5);
plot(phx,ones(length(phx))*ph,'--k','LineWidth',1.5)
set(gca,'Xlim',[0 3])
title('Cross-Shore')
xlabel('V/V_r_m_s')
ylabel('Height Above Bed (cm)')

subplot(133)
Wnorm = W./Wrms;
Wbar = mean(Wnorm,1);
stdv = std(Wnorm,0,1);
p(1) = plot(Wbar,rbcm,...
    'Color','k','Marker','d','MarkerFaceColor','k','LineWidth',1.5);hold on
p(2) = plot((Wbar-stdv),rbcm,...
    'Color',[0.6 0.6 0.6],'LineStyle','--','LineWidth',1.5);
p(3) = plot((Wbar+stdv),rbcm,...
    'Color',[0.6 0.6 0.6],'LineStyle','--','LineWidth',1.5);
plot(phx,ones(length(phx))*ph,'--k','LineWidth',1.5)
set(gca,'Xlim',[-0.5 0.5])
title('Vertical')
xlabel('W/W_r_m_s')
ylabel('Height Above Bed (cm)')
suptitle('AD5117 - Above Pneumatophores')

%try plotting RMS normalized velocities by z/hc (measurement height/height
%of canopy)
f10 = figure(10);
set(f10,'PaperOrientation','portrait',...
    'position',[200 200   800   600]);

phx = linspace(-5,5,10);
plot(phx,ones(length(phx)),'--','Color',[0 0 0],'LineWidth',1.5), hold on
zhc = dat.a6.rbcm/ph;
Unorm = dat.a6.U./dat.a6.Urms;
Ubar = mean(Unorm,1); 
p(1) = plot(Ubar,zhc,...
    'Color','k','Marker','x','MarkerFaceColor',[0.5 0.5 0.5],'Color',[0.5 0.5 0.5],'MarkerSize',12,...
    'LineWidth',1.5);
Vnorm = dat.a6.V./dat.a6.Vrms;
Vbar = mean(Vnorm,1);
p(2) = plot(Vbar,zhc,...
    'Color','k','Marker','d','MarkerFaceColor',[0.5 0.5 0.5],'Color',[0.5 0.5 0.5],'LineWidth',1.5);
Wnorm = dat.a6.W./dat.a6.Wrms;
Wbar = mean(Wnorm,1);
p(3) = plot(Wbar,zhc,...
    'Color','k','Marker','o','MarkerFaceColor',[0.5 0.5 0.5],'Color',[0.5 0.5 0.5],'LineWidth',1.5);
zhc = dat.a7.rbcm/ph;
Unorm = dat.a7.U./dat.a7.Urms;
Ubar = mean(Unorm,1); 
p(4) = plot(Ubar,zhc,...
    'Color','k','Marker','x','MarkerFaceColor',[0 0 0],'Color',[0 0 0],'MarkerSize',12,...
    'LineWidth',1.5);
Vnorm = dat.a7.V./dat.a7.Vrms;
Vbar = mean(Vnorm,1);
p(5) = plot(Vbar,zhc,...
    'Color','k','Marker','d','MarkerFaceColor',[0 0 0],'Color',[0 0 0],'LineWidth',1.5);
Wnorm = dat.a7.W./dat.a7.Wrms;
Wbar = mean(Wnorm,1);
p(6) = plot(Wbar,zhc,...
    'Color','k','Marker','o','MarkerFaceColor',[0 0 0],'Color',[0 0 0],'LineWidth',1.5);
set(gca,'Xlim',[-1 1.5],'Ylim',[0.4 1.2],'FontSize',12,'FontName','Arial')
xlabel('U/U_r_m_s')
ylabel('z/h_c')
leg = legend(p(4:6),{'Along-Shore','Cross-Shore','Vertical'});
set(leg,'box','off','position',[0.21 0.83 0.05 0.05])

