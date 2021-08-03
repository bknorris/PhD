%Calculate Uc from Longuet-Higgens (1953) to check if the wave climate of
%the mudflat-mangrove system is enough  to induce drift through wave
%stress.

clear
load('c:\Users\bkn5\Projects\Mekong_W2015\DataAnalysis\Environmental\V5108_HTA.mat')
load('c:\Users\bkn5\Projects\Mekong_W2015\DataAnalysis\Paper1\HTAtimes.mat')
zp = 0.196; %pressure sensor height above bottom
zuv = 0.380; %height of velocity sensor above bottom
intv = 2; %averaging interval (minutes)
heading = 110; %rotate to x-shore azimuth
start = HTA.times.t1;
stop = HTA.times.e1;
savedatfile = 1;
datdir = 'c:\Users\bkn5\Projects\Mekong_W2015\DataAnalysis\Spectra\Paper1\';

%crop data fields to the specified start/stop times of the experiment
ind = find(ADV.datetime >= start & ADV.datetime <=stop);
time = ADV.datetime(ind);
pressure = ADV.Pres(ind);
u = ADV.U(ind);v = ADV.V(ind);
fs = ADV.Metadata.instmeta.samprate; %sampling frequency for ADVs
x = (sin(ADV.Metadata.inst_lat/57.29578))^2;

%Iterate through time intervals
avt = fs*intv*60; %# samples in interval
idx = avt:avt:length(u);
idxx = [1 idx];    
it = 1;
while it < length(idxx)
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
    U(isnan(U)) = 0;U = detrend(U);
    V(isnan(V)) = 0;V = detrend(V);
    [m,n] = size(U);
    
    %rotate to the cross-shore direction
    rot = heading*pi/180;
    U = U.*(ones(size(U))*cos(rot)) + ...
        V.*(ones(size(V))*sin(rot));
    V = -U.*(ones(size(U))*sin(rot)) + ...
        V.*(ones(size(V))*cos(rot));
    
    %Spectra settings
    nfft = 0.25*m;
    window = hanning(m,'periodic');
    noverlap = floor(0.7*length(window));
    min_f = 0.01; % mininum frequency, below which no correction is applied
    max_f = 1.5; % maximum frequency, above which no correction is applied
    
    g = 9.81;
    rho = 1.009629812357770e+03;
    for i = 1:length(p)
        g(i,:) = 9.780318*(1.0+(5.2788E-3+2.36E-5*x)*x)+1.092E-6*p(i,:);
        h(i,:) = ((((-1.82E-15*p(i,:)+2.279E-10)*p(i,:)-2.2512E-5)*p(i,:)+9.72659)*p(i,:))/g(i,:);
    end
    h = mean(h);

    % Determine velocity spectra for u and v
    [Guu, f] = pwelch(U,window,noverlap,nfft,fs);
    [Gvv, ~] = pwelch(V,window,noverlap,nfft,fs);
    df = f(3)-f(2);
    
    % determine wave number
    omega = 2*pi.*f;
    k = qkhf(omega,h)./h;
    
    %compute linear wave transfer function
    kh = k*h;
    nf = length(omega);
    Huv = ones(nf,1);
    
    % combine horizontal velocity spectra
    Guv = Guu + Gvv;
    
    % Determine wave height for velocity spectra
    Snu = Guv./(Huv.^2);
    
    % create cut off freqency, so noise is not magnified
    ff = min(find(f>=min_f));
    lf = max(find(f<=max_f));
    
    % Determine rms wave height (mult by another sqrt(2) for Hs)
    % Thornton and Guza say Hrms = sqrt(8 mo)
    Hrmsu = 2*sqrt( 2*sum( Snu(ff:lf)*df ) );
    Hs = Hrmsu*sqrt(2);
    
    %calculate L-H parameters: A, B, C, D
    z = zuv;
    a = Hs/2;
    
%     Uc = (3/4).*(k*a).*((a.*omega)./(sinh(kh.^2)));
    Uu = (5*(a^2).*omega.*k)./(4*sinh(kh.^2));
    Uc(it,:) = max(Uu);
    T(it,:) = ts;
    it = it+1;
end

%load the VecPros, compare mean bottom currents to Uc
load('C:\Users\bkn5\Projects\Mekong_W2015\DataAnalysis\Paper1\HTAday1Vels.mat')
fs = 50;
intv = 2;
heading = 10;
avt = fs*intv*60; %# samples in interval
idx = avt:avt:length(HTA.vpro3.x);
idxx = [1 idx];    
it = 1;
while it < length(idxx)
    indx = idxx(it):idxx(it+1);
    x = HTA.vpro3.x(indx,:);
    y = HTA.vpro3.y(indx,:);
    v = x(:,3); 
    u = y(:,3);
    
    %rotate to the cross-shore direction
    rot = heading*pi/180;
    u = u.*(ones(size(u))*cos(rot)) + ...
        v.*(ones(size(v))*sin(rot));
    v = -u.*(ones(size(u))*sin(rot)) + ...
        v.*(ones(size(v))*cos(rot));
    [spd,dir] = cmguv2spd(x(:,3),y(:,3));
    S(it,:) = mean(spd);
    Dir(it,:) = mean(dir);
    Uvp(it,:) = nanmean(u);
    it = it+1;
end
if length(Uvp) ~= length(Uc)
    b = length(Uvp);c = length(Uc);
    Uvp(b:c) = 0;
end

figure;
q(1) = plot(T,Uc,'k','LineWidth',1.5,'Marker','o','MarkerFaceColor',[1 1 1]);
hold on
q(2) = plot(T,Uvp.*100,'r','LineWidth',1.5,'Marker','^','MarkerFaceColor',[1 1 1]);
datetick('x','HH:MM','keepticks','keeplimits')
ylabel('Current Velocity (cm/s)','FontSize',14,'FontName','Cambria')
xlabel('Time on 07/03/2015','FontSize',14,'FontName','Cambria')
leg = legend(q,{'Uc (L-H 1953)';'Mean Current (x = 20 cm)'});set(leg,'box','off')
set(gca,'LineWidth',1.5,'FontSize',14,'FontName','Cambria')

disp(['Uc (L-H 1953): ' num2str(mean(Uc)) ' cm/s'])
disp(['Mean Current measured at x = 20cm: ' num2str(mean(Uvp*100)) ' cm/s'])