clear
load('C:\Users\bkn5\Projects\Mekong_W2015\Data\Aquadopp\FSS\AD5116_12March2015.mat')
name = 'HTAday3';
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
dat.a6.rbcm = rbcm;
sr = 8;
intv = 2; %minutes
avt = intv*sr*60;
ind = [1 avt:avt:length(w)];

%compute speed and direction over the interval
[~,n] = size(u);m = length(ind)-1;
Spd = zeros(m,n);Dir = zeros(m,n);T = zeros(m,1);
for i = 1:m
    for ii = 1:n
        U = u(ind(i):ind(i+1),ii);
        V = v(ind(i):ind(i+1),ii);
        [dat.a6.Spd(i,ii),dat.a6.Dir(i,ii)] = cmguv2spd(nanmean(U),nanmean(V));
        dat.a6.T(i,:) = t(ind(i));
    end
end

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
dat.a7.rbcm = rbcm;
sr = 8;
avt = intv*sr*60;
ind = [1 avt:avt:length(w)];

%compute speed and direction over the interval
[~,n] = size(u);m = length(ind)-1;
Spd = zeros(m,n);Dir = zeros(m,n);T = zeros(m,1);
for i = 1:m
    for ii = 1:n
        U = u(ind(i):ind(i+1),ii);
        V = v(ind(i):ind(i+1),ii);
        [dat.a7.Spd(i,ii),dat.a7.Dir(i,ii)] = cmguv2spd(nanmean(U),nanmean(V));
        dat.a7.T(i,:) = t(ind(i));
    end
end

%%%%%%%% Plot routine

%plot time-averaged depth profiles of speed, direction
f1 = figure(1);
set(f1,'PaperOrientation','portrait',...
    'position',[200 200   1000   600]);
set(gcf,'color','w','PaperPositionMode','auto')
pq(1) = subplot(121);
%pneumatophore height of canopy
phx = linspace(-20,20,10);
plot(phx,ones(length(phx))*ph,'--k','LineWidth',1.5), hold on

%bare bed
mSpd = mean(dat.a6.Spd);Spdstd = std(dat.a6.Spd);

h1 = herrorbar(mSpd*100,dat.a6.rbcm,2*100*Spdstd);
set(h1,'LineWidth',1.5,'Color',[0.5 0.5 0.5]), hold on
plot(mSpd*100,dat.a6.rbcm,'-s','Color',[0.5 0.5 0.5],'LineWidth',2,...
    'markersize',10,'markerfacecolor',[1 1 1])

%above canopy
mSpd = mean(dat.a7.Spd);Spdstd = std(dat.a7.Spd);

h2 = herrorbar(mSpd*100,dat.a7.rbcm,2*100*Spdstd);
set(h2,'LineWidth',1.5,'Color','k')
plot(mSpd*100,dat.a7.rbcm,'-^','Color','k','LineWidth',2,...
    'markersize',10,'markerfacecolor',[1 1 1])
set(gca,'YDir','normal','XLim',[-5 15],...
    'FontSize',14,'FontName','Cambria','LineWidth',1.5)
xlabel('Speed (cm/s)','FontSize',14,'FontName','Cambria')
ylabel('Height Above Bed (cm)','FontSize',14,'FontName','Cambria')

pq(2) = subplot(122);
%pneumatophore height of canopy
phx = linspace(0,360,10);
plot(phx,ones(length(phx))*ph,'--k','LineWidth',1.5), hold on

%bare bed
mDir = mean(dat.a6.Dir);
plot(mDir,dat.a6.rbcm,'-s','Color',[0.5 0.5 0.5],'LineWidth',2,...
    'markersize',10,'markerfacecolor',[1 1 1])

%above canopy
mDir = mean(dat.a7.Dir);Dirstd =std(dat.a7.Dir);
plot(mDir,dat.a7.rbcm,'-^','Color','k','LineWidth',2,...
    'markersize',10,'markerfacecolor',[1 1 1])
set(gca,'YDir','normal','XLim',[0 180],'XTick',0:45:180,...
    'YTickLabel',[],'FontSize',14,'FontName','Cambria','LineWidth',1.5)
xlabel('Direction (Deg)','FontSize',14,'FontName','Cambria')
hold off

set(pq(1),'position',[0.1 0.1 0.4 0.85])
set(pq(2),'position',[0.54 0.1 0.4 0.85])

export_fig(['c:\Users\bkn5\Projects\Mekong_W2015\Figures\Paper1\' 'AQDPsSpdDir' name],'-jpg')