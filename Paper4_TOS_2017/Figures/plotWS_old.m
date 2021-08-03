%calculate velocity and wave direction from an offshore ADV - Overview for
%complete Time Series. Load V5108 and V5109 for the HTA and VTA,
%respectively

clear %a good idea
close all
load('c:\users\bkn5\Projects\Mekong_W2015\Data\Vector\SW_Mudflat\V5108_030315.mat')
name = 'V5108'; %instrument
plotfigs = 0;
zp = 0.196; %pressure sensor height above bottom
zuv = 0.380; %height of velocity sensor above bottom
intv = 1; %averaging interval (minutes)
start = datenum(2015,03,06,12,00,00);
stop = datenum(2015,03,12,09,00,00);
% savedatfile = 1;
% figdir = 'c:\Users\bkn5\Projects\Mekong_W2015\Figures\Spectra&Waves\';
% datdir = 'c:\Users\bkn5\Projects\Mekong_W2015\DataAnalysis\Spectra\FSS3\';

%crop data fields to the specified start/stop times of the experiment
ind = find(ADV.datetime >= start & ADV.datetime <=stop);
pressure = ADV.Pres(ind);
u = ADV.U(ind);v = ADV.V(ind);
fs = ADV.Metadata.instmeta.samprate; %sampling frequency for ADVs

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
    U = u(ind);V = v(ind);
    if rem(length(p), 2) == 1
        n = length(p)-1; %force inputs to be even
        p = p(1:end-1);
        U = U(1:end-1);
        V = V(1:end-1);
    else
        n = length(p);
    end
    nfft = n/2;
    window = hanning(n/5);
    min_f = 0.05; % mininum frequency, below which no correction is applied
    max_f = 0.5; % maximum frequency, above which no correction is applied
    
    nlin = 1E2;
    maxgaps = 1E5;
    p = cmgbridge(p,nlin,maxgaps,maxgaps);
    U = cmgbridge(U,nlin,maxgaps,maxgaps);
    V = cmgbridge(V,nlin,maxgaps,maxgaps);
    rho = 1.009629812357770e+03;
    
    ws = puv(p,U,V,h,zp,zuv,fs,nfft,rho,window,min_f,max_f);
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

    title(['\bf\it' datestr(start,'dd/mm/yyyy')],'FontSize',14) 
    s(1) = subplot(311);
    p(1) = plot(t,wvstats.depth,'-','Color','k','LineWidth',1.5);hold on
    ylabel('\bf\itDepth (m)','FontSize',14)
    set(gca,'XTickLabel',[],'LineWidth',1.5,'FontSize',14,'XLim',[0 length(t)])
    
    s(2) = subplot(312);
    p(2) = plot(t,wvstats.hrmsp,'-^','Color','k','LineWidth',1.5);hold on
    p(3) = plot(t,wvstats.hrmsuv,'-o','Color',[0.65 0.65 0.65],'LineWidth',1.5);
    ylabel('\bf\itHrms (m)','FontSize',14)
    set(gca,'XTickLabel',[],'LineWidth',1.5,'FontSize',14,'XLim',[0 length(t)])
    
    l(1) = legend('From Pressure','From Velocity');
    
    s(3) = subplot(313);
    p(4) = plot(t,wvstats.udir,'d','Color','k','MarkerSize',8,'MarkerFaceColor','k');hold on
    p(5) = plot(t,(360-wvstats.azm),'o','Color',[0.65 0.65 0.65],'MarkerSize',8,'MarkerFaceColor',[0.65 0.65 0.65]);
    ylabel('\bf\itF_m (Wm^-^1)','FontSize',14)
    set(gca,'LineWidth',1.5,'FontSize',14,'XLim',[0 length(t)],'YLim',[0 400],'YTick',0:90:360)
    xlabel('\bf\itTime Elapsed (min)','FontSize',14)
    ylabel('\bf\itGeographic Direction (deg)','FontSize',14)
    l(2) = legend('Current direction','Wave direction');
    
    set(s(1),'position',[0.1 0.66 0.82 0.28])
    set(s(2),'position',[0.1 0.38 0.82 0.28])
    set(s(3),'position',[0.1 0.1 0.82 0.28])
    
    fpath = figdir;fname = [dname '_wavestats'];
    export_fig([fpath fname],'-png','-native')
    disp(['Figure ' fname '.png saved'])
end
