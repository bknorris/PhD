%calculate velocity and wave direction from an offshore Aquadopp

clear %a good idea
close all
load('d:\Projects\Mekong_F2014\Data\Aquadopp\F2F\HR4.mat')
name = 'HR4'; %instrument
dname = 'F2F';
zp = 0.082; %pressure sensor height above bottom (m)
zuv = 0.082; %height of velocity sensor above bottom (m)
intv = 600; %averaging interval (seconds)
start = datenum(2014,09,29,14,50,00);
stop = datenum(2014,09,30,08,30,00);
savedatfile = 1;
datdir = 'd:\Projects\Mekong_F2014\DataAnalysis\Paper2\Spectra\';

%crop data fields to the specified start/stop times of the experiment
ind = find(aqdp.datenum >= start & aqdp.datenum <=stop);
time = aqdp.datenum(ind);
pressure = aqdp.pressure(ind);
u = nanmean(aqdp.u(ind),2);v = nanmean(aqdp.v(ind),2); %depth average profile
fs = str2double(aqdp.metadata.samprate(1)); %sampling frequency for aqdps
x = (sin(aqdp.metadata.lat/57.29578))^2;

%prep data for analysis by spectrally filling NaNs with cmgbridge
u = cmgbridge(u,100,100,1000);
v = cmgbridge(v,100,100,1000);
pressure = cmgbridge(pressure,100,100,1000);

%Set up wavestat structure
wvstats = struct();
wvstats.info.start = [datestr(start,'dd-mm-yyyy HH:MM:SS') ' ICT'];
wvstats.info.stop = [datestr(stop,'dd-mm-yyyy HH:MM:SS') ' ICT'];
wvstats.info.depname = dname;
wvstats.info.inst = name;
wvstats.info.intv = [num2str(intv) ' minutes'];
wvstats.info.fieldinfo = {'Structure wvstats field descriptions:';...
   'wvstats.uspd = current velocity';...
   'wvstats.udir = current direction';...
   'wvstats.u = [mean] easterly current velocity';...
   'wvstats.v = [mean] northerly current velocity';...
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

%Iterate through time intervals
avt = fs*intv; %# samples in interval
idx = avt:avt:length(u);
idxx = [1 idx];    
it = 1;
while it < length(idxx)
    if (it/1000)-floor(it/1000)==0
        disp(['Iteration number ' num2str(it)])
    end
    indx = idxx(it):idxx(it+1);
    p = pressure(indx);
    U = u(indx);
    V = v(indx);
    ts = time(idxx(it)); %save timestamp of the spectra interval
    if rem(length(p), 2) == 1
        n = length(p)-1; %force inputs to be even
        p = p(1:end-1);
        U = U(1:end-1);
        V = V(1:end-1);
    else
        n = length(p);
    end
    U(isnan(U)) = 0; %deal with leading/trailing NaNs
    V(isnan(V)) = 0;
    p(isnan(p)) = 0;
    
    %Spectra settings
    nfft = n/5;
    window = hanning(n/5);
    min_f = 0.05; % mininum frequency, below which no correction is applied
    max_f = 1; % maximum frequency, above which no correction is applied
    
    g = 9.81;
    rho = 1.009629812357770e+03;
    for i = 1:length(p)
        g(i,:) = 9.780318*(1.0+(5.2788E-3+2.36E-5*x)*x)+1.092E-6*p(i,:);
        h(i,:) = ((((-1.82E-15*p(i,:)+2.279E-10)*p(i,:)-2.2512E-5)*p(i,:)+9.72659)*p(i,:))/g(i,:);
    end
    h = mean(h);
    
    ws = puv(p,U,V,h,zp,zuv,fs,nfft,rho,window,min_f,max_f);
    [spd,dir] = cmguv2spd(U,V);
    %if there are more than half zeros in the velocity (because the
    %instrument is out of the water), set the azimuth and direction to NaN
    %to exclude it from the dataset.
    
    if sum(U == 0) > length(U)/2 || sum(V == 0) > length(V)/2
        ws.azr = NaN;
        dir = NaN;
    end
        
    wvstats.uspd(it,1) = mean(spd);
    wvstats.udir(it,1) = mean(dir);
    wvstats.u(it,1) = nanmean(U);
    wvstats.v(it,1) = nanmean(V);
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
    wvstats.time(it,1) = ts;
    
    it = it+1;
end
wvstats.lfcutoff = min_f;
wvstats.hfcutoff = max_f;

if savedatfile
    %save file
    filename = [name '_' dname];
    sfname = ['WaveStats' filename];
    save([datdir sfname],'wvstats','-v7.3')
    disp(['Wave Statistics data file saved as ' sfname])
end
% % if plotfigs
%     f1 = figure(1);
%     set(f1,'PaperOrientation','portrait',...
%     'position',[400 200   1000   800]);
%     set(gcf,'color','w','PaperPositionMode','auto')
% 
%     title(['\bf\it' datestr(start,'dd/mm/yyyy')],'FontSize',14) 
%     s(1) = subplot(311);
%     p(1) = plot(t,wvstats.depth,'-','Color','k','LineWidth',1.5);hold on
%     ylabel('\bf\itDepth (m)','FontSize',14)
%     set(gca,'XTickLabel',[],'LineWidth',1.5,'FontSize',14,'XLim',[0 length(t)])
%     
%     s(2) = subplot(312);
%     p(2) = plot(t,wvstats.hrmsp,'-^','Color','k','LineWidth',1.5);hold on
%     p(3) = plot(t,wvstats.hrmsuv,'-o','Color',[0.65 0.65 0.65],'LineWidth',1.5);
%     ylabel('\bf\itHrms (m)','FontSize',14)
%     set(gca,'XTickLabel',[],'LineWidth',1.5,'FontSize',14,'XLim',[0 length(t)])
%     
%     l(1) = legend('From Pressure','From Velocity');
%     
%     s(3) = subplot(313);
%     p(4) = plot(t,wvstats.udir,'d','Color','k','MarkerSize',8,'MarkerFaceColor','k');hold on
%     p(5) = plot(t,(360-wvstats.azm),'+','Color',[0.65 0.65 0.65],'MarkerSize',8,'MarkerFaceColor',[0.65 0.65 0.65]);
%     ylabel('\bf\itF_m (Wm^-^1)','FontSize',14)
%     set(gca,'LineWidth',1.5,'FontSize',14,'XLim',[0 length(t)],'YLim',[0 400],'YTick',0:90:360)
%     xlabel('\bf\itTime Elapsed (min)','FontSize',14)
%     ylabel('\bf\itGeographic Direction (deg)','FontSize',14)
%     l(2) = legend('Current direction','Wave direction');
%     
%     set(s(1),'position',[0.1 0.66 0.82 0.28])
%     set(s(2),'position',[0.1 0.38 0.82 0.28])
%     set(s(3),'position',[0.1 0.1 0.82 0.28])
%     
%     fpath = figdir;fname = [dname '_wavestats'];
%     export_fig([fpath fname],'-png','-native')
%     disp(['Figure ' fname '.png saved'])
% % end
