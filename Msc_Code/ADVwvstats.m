%calculate velocity and wave direction from an offshore ADV

clear %a good idea
close all
load('c:\users\bkn5\Projects\Mekong_W2015\Data\Vector\DPS2\V5109_130315.mat')
dname = 'DPS2'; %deployment
name = 'V5109'; %instrument
plotfigs = 1;
zp = 0.120; %pressure sensor height above bottom
zuv = 0.290; %height of velocity sensor above bottom
intv = 1; %averaging interval (minutes)
% start = datenum(2015,03,07,14,20,00);
% stop = datenum(2015,03,07,16,00,00);
% start = datenum(2015,03,08,15,00,00);
% stop = datenum(2015,03,08,18,45,00);
% start = datenum(2015,03,10,15,31,00);
% stop = datenum(2015,03,10,16,52,00);
start = datenum(2015,03,14,06,40,00);
stop = datenum(2015,03,14,09,08,00);
savedatfile = 1;
figdir = 'c:\Users\bkn5\Projects\Mekong_W2015\Figures\Spectra&Waves\';
datdir = 'c:\Users\bkn5\Projects\Mekong_W2015\DataAnalysis\Spectra\DPS2\';

%crop data fields to the specified start/stop times of the experiment
ind = find(ADV.datetime >= start & ADV.datetime <=stop);
pressure = ADV.Pres(ind);
u = ADV.U(ind);v = ADV.V(ind);
fs = ADV.Metadata.instmeta.samprate; %sampling frequency for ADVs
ind2 = find(ADV.Sensor.Datetime >= start & ADV.Sensor.Datetime <=stop);
temp = ADV.Sensor.Temp(ind2); %Vector records temperature at 1Hz
temperature = spline(ADV.Sensor.Datetime(ind2),temp,ADV.datetime(ind));
x = (sin(ADV.Metadata.inst_lat/57.29578))^2;

%loop over the specified sample interval
avt = fs*intv*60; %# samples in interval
idx = avt:avt:length(u);
idxx = [1 idx];
wvstats = struct();
wvstats.info.start = [datestr(start,'dd-mm-yyyy HH:MM:SS') ' ICT'];
wvstats.info.stop = [datestr(stop,'dd-mm-yyyy HH:MM:SS') ' ICT'];
wvstats.info.depname = dname;
wvstats.info.inst = name;
wvstats.info.intv = [num2str(intv) ' minutes'];
wvstats.info.fieldinfo = {'Structure wvstats field descriptions:';...
   'wvstats.uspd = current velocity';...
   'wvstats.udir = current direction';...
   'wvstats.hrmsp = Hrms (=Hmo) from pressure';...
   'wvstats.hrmsuv = Hrms from u,v';...
   'wvstats.rorb = Representative orbital velocity amplitude in freq. band ( lf_cutoff <= f <= hf_cutoff ) (m/s)';...
   'wvstats.omegar = Representative orbital velocity (radian frequency)';...
   'wvstats.Tr = Representative orbital velocity period (s)';...
   'wvstats.phir = Representative orbital velocity direction (angles from x-axis, positive ccw)';...
   'wvstats.azm = Represntative orb. velocity direction (deg; geographic azimuth; [0 - 360], 0 = North; ambiguous around 180 degrees)';...
   'wvstats.orblo = orb in freq. band (f <= lf_cutoff) (m/s)';...
   'wvstats.orbhi = orb in freq. band (f >= hf_cutoff) (m/s)';...
   'wvstats.orbifg = orb in infragravity freq. band (lf_cutoff f <= 1/20) (m/s)'};

    
it = 1;
%Iterate through time intervals %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
while it < length(idxx)
    disp(['Iteration number ' num2str(it)])
    ind = idxx(it):idxx(it+1);
    p = pressure(ind);
    t = temperature(ind);
    U = u(ind);V = v(ind);
    if rem(length(p), 2) == 1
        n = length(p)-1; %force inputs to be even
        p = p(1:end-1);
        t = t(1:end-1);
        U = U(1:end-1);
        V = V(1:end-1);
    else
        n = length(p);
    end
    nfft = n/2;
    min_f = 0.05; % mininum frequency, below which no correction is applied
    max_f = 0.5; % maximum frequency, above which no correction is applied
    
    nlin = 1E2;
    maxgaps = 1E5;
    p = cmgbridge(p,nlin,maxgaps,maxgaps);
    t = cmgbridge(t,nlin,maxgaps,maxgaps);
    U = cmgbridge(U,nlin,maxgaps,maxgaps);
    V = cmgbridge(V,nlin,maxgaps,maxgaps);
    
    salt = repmat(20,length(t),1); %assume constant salinity of 20PSU
    rho = SeaDensity(salt,t,p);
%     Calculate Depth (h) from the pressure signal
%     From: UNESCO (1983): Algorithms for computation of fundamental properties
%     of seawater. UNESCO technical papers in marine science 44:1-55.
    g = zeros(length(p),1);h = zeros(length(p),1);
    for i = 1:length(p)
        g(i,:) = 9.780318*(1.0+(5.2788E-3+2.36E-5*x)*x)+1.092E-6*p(i,:);
        h(i,:) = ((((-1.82E-15*p(i,:)+2.279E-10)*p(i,:)-2.2512E-5)*p(i,:)+9.72659)*p(i,:))/g(i,:);
    end
    rho = mean(rho);
    h = mean(h);
    
    ws = puv(p,U,V,h,zp,zuv,fs,nfft,rho,hanning(n),min_f,max_f);
    [spd,dir] = cmguv2spd(U,V);
    wvstats.uspd(it,1) = mean(spd);
    wvstats.udir(it,1) = mean(dir);
    wvstats.hrmsp(it,1) = ws.Hrmsp;
    wvstats.hrmsuv(it,1) = ws.Hrmsu;
    wvstats.rorb(it,1) = ws.ubr;
    wvstats.omegar(it,1) = ws.omegar;
    wvstats.Tr(it,1) = ws.Tr;
    wvstats.phir(it,1) = ws.phir;
    wvstats.azm(it,1) = ws.azr;
    wvstats.orblo(it,1) = ws.ublo;
    wvstats.orbhi(it,1) = ws.ubhi;
    wvstats.orbifg(it,1) = ws.ubig;
    wvstats.depth(it,1) = mean(p)+zp;
    
    it = it+1;
end
wvstats.lfcutoff = min_f;
wvstats.hfcutoff = max_f;
t = 1:length(wvstats.uspd);
transect = 280*ones(length(t),1);
if savedatfile
    %save file
    filename = [name '_' dname];
    sfname = ['WaveStats' filename];
    save([datdir sfname],'wvstats','-v7.3')
    disp(['Drag force data file saved as ' sfname])
end
if plotfigs
    f1 = figure(1);
    set(f1,'PaperOrientation','portrait',...
    'position',[400 200   1000   800]);
    set(gcf,'color','w','PaperPositionMode','auto')
    
    s(1) = subplot(311);
    p(1) = plot(t,wvstats.depth,'-','Color','k','LineWidth',1.5);hold on
    ylabel('\bf\itDepth (m)','FontSize',14)
    set(gca,'XTickLabel',[],'LineWidth',1.5,'FontSize',14,'XLim',[0 length(t)],'YLim',[0 2],'YTick',0:0.5:2)
    title(['\bf\it' datestr(start,'dd/mm/yyyy')],'FontSize',14) 
    
    s(2) = subplot(312);
    p(2) = plot(t,wvstats.hrmsp,'-^','Color','k','LineWidth',1.5);hold on
    p(3) = plot(t,wvstats.hrmsuv,'-o','Color',[0.65 0.65 0.65],'LineWidth',1.5);
    ylabel('\bf\itHrms (m)','FontSize',14)
    set(gca,'XTickLabel',[],'LineWidth',1.5,'FontSize',14,'XLim',[0 length(t)],'YLim',[0 0.6],'YTick',0:0.25:0.5)
    l(1) = legend('From Pressure','From Velocity');
    set(l(1),'LineWidth',1.5,'position',[0.16 0.54 0.1 0.1])
    
    s(3) = subplot(313);
    p(4) = plot(t,wvstats.udir,'d','Color','k','MarkerSize',8,'MarkerFaceColor','k');hold on
    p(5) = plot(t,(360-wvstats.azm),'o','Color',[0.65 0.65 0.65],'MarkerSize',8,'MarkerFaceColor',[0.65 0.65 0.65]);
    p(6) = line(t,transect,'Color','r','LineWidth',1.5);
    ylabel('\bf\itF_m (Wm^-^1)','FontSize',14)
    set(gca,'LineWidth',1.5,'FontSize',14,'XLim',[0 length(t)],'YLim',[90 400],'YTick',90:90:360)
    xlabel('\bf\itMinutes After Low Tide','FontSize',14)
    ylabel('\bf\itAzimuth (deg)','FontSize',14)
    l(2) = legend('Current direction','Wave direction');
    set(l(2),'LineWidth',1.5,'position',[0.16 0.12 0.1 0.1])
    
    set(s(1),'position',[0.1 0.66 0.82 0.28])
    set(s(2),'position',[0.1 0.38 0.82 0.28])
    set(s(3),'position',[0.1 0.1 0.82 0.28])
    
    
    fpath = figdir;fname = [name '_' dname '_wavestats'];
    export_fig([fpath fname],'-png','-native')
    disp(['Figure ' fname '.png saved'])
end
