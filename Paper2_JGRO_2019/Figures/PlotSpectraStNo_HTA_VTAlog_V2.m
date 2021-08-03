%Plot cross-shore velocity spectra from the 3 VPs of the HTA and the 3 VPs
%of the VTA experiment. Input time series are 30 minutes of data during the
%individual experiments.
clear

figdir = 'f:\GradSchool\DataAnalysis\Paper2\WorkingFigures\Spectra\';
datdir = 'd:\Mekong_W2015\DataAnalysis\Paper1\';
files = {'HTA_1Vels.mat';'HTA_2Vels.mat';'VTA_2vp1Vels.mat';};
days = {'day1';'day2';'day4'};
txt = {'x = 20 cm';'x = 10 cm';''};
%%%Plot Routine

f1 = figure(1);
set(f1,'PaperOrientation','portrait',...
    'position',[400 100   1000   400],...
    'renderer','painters');
line = {'-';'--';'-.'};
heading = [20 20 96];
peaks = [2.3 2.4 1.9];
vps = [3 2 1];
sp = zeros(3,1);
data = struct();
for i = 1:3
    load([datdir files{i}])
    fn = fieldnames(dat);
    if i == 2
        bins = 9:23;
    else
        bins = 1:5;
    end
    %Spectra settings
    fs = 50;
    win = 600; %seconds (10 minutes)
    nwin = fs*win;
    swin = fs*40; %20 second averaging window (to smooth)
    DOF = round((nwin/swin)*2);
    if i == 1
        disp(['50% Hamming windowed spectra with ' sprintf('%0.0f',DOF) ' degrees of freedom'])
    end
    
    ii = vps(i);
    sp(i) = subplot(1,3,i);
    x = mean(dat.(fn{ii}).x(:,bins),2);
    y = mean(dat.(fn{ii}).y(:,bins),2);
    rot = (pi*heading(i))/180;
    T = [cos(rot) -sin(rot);...
        sin(rot) cos(rot)];
    vels = [x y];
    V = vels*T';
    x = V(:,1);y = V(:,2);
    %Calculate difference between a single bin in the middle of the profile
    mid = median(bins);
    xm1 = dat.(fn{ii}).x(:,mid+1);
    xm2 = dat.(fn{ii}).x(:,mid);
    xdif = xm1-xm2;
    
    %calculate RMS velocities (Luhar et al. 2013)
    Ec = (1/length(y))*sum(y);Nc = (1/length(x))*sum(x);
    Ewrms = sqrt((1/length(y))*sum((y-Ec).^2));
    Nwrms = sqrt((1/length(x))*sum((x-Nc).^2));
    Uc = sqrt(Ec^2+Nc^2);Uwrms = sqrt(Ewrms^2+Nwrms^2);
    Uw = sqrt(2)*Uwrms;
    %calculate Eulerian time scale (autocorr of x) %%Not sure if this is
    %correct. Check Tanino & Nepf (2008) if it is necessary to compute for
    %other deployments.
    [acor,lags] = xcov(x,'coeff');
    acor = acor(length(x):end);
    id = find(acor <= 0,1,'first');
    tu = trapz(acor(1:id-1))/fs;
    %save data to structure
    data.(days{i}).(fn{ii}).Uc = Uc;
    data.(days{i}).(fn{ii}).Uw = Uw;
    data.(days{i}).(fn{ii}).tu = tu;
    %calculate psd
    x = x(1:nwin*6);
    xdif = xdif(1:nwin*6);
    [Cuu,F,cf] = pwelch(detrend(x),hanning(swin),swin*0.5,nwin,fs,'confidencelevel',0.95);
    [Cdif,~,cfdif] = pwelch(detrend(xdif),hanning(swin),swin*0.5,nwin,fs,'confidencelevel',0.95);
    Cuu(F < 0.05) = NaN;
    Cdif(F < 0.05) = NaN;
    cfdif(F < 0.05,:) = [];
    cf(F < 0.05,:) = [];
    ff = F;
    ff(F < 0.05,:) = [];
    %plot 5/3 slope on figures
    if i == 3
        int = 0.0015;
        text(3,5E-4,'^-^5^/^3'),hold on
    else
        int = 0.005;
        text(3,1E-3,'^-^5^/^3'),hold on
    end
    xs = linspace(2,45,length(Cuu));
    ys = int.*(xs.^(-5/3));
    plot(xs,ys,'Color','r','LineWidth',1.5);hold on
    xx = [ff;flipud(ff)];    
    yy = [cfdif(:,1);flipud(cfdif(:,2))];
    fill(xx,yy,[0.7 0.7 0.7],...
        'EdgeColor','none'), hold on
    plot(F,Cdif,'--',...
        'Color',[0.4 0.4 0.4],...
        'LineWidth',1.5)
    yy = [cf(:,1);flipud(cf(:,2))];
    fill(xx,yy,[0.7 0.7 0.7],...
        'EdgeColor','none'), hold on
    plot(F,Cuu,...
        'Color','k',...
        'LineWidth',1.5)
    pid = find(F >= peaks(i),1,'first');
    text(2,4E-2,txt{i},'FontSize',14)
    set(gca,'yscale','log','xscale','log')
end

%plot adjustments
set(sp,'xlim',[0.1 10])
set([sp(1) sp(2)],'ylim',[1E-5 1E-1])
set(sp(3),'ylim',[1E-6 1E-2])
set(sp(1),'position',[0.1 0.15 0.26 0.7])
set(sp(2),'position',[0.4 0.15 0.26 0.7],...
    'yticklabel',[])
set(sp(3),'position',[0.72 0.15 0.26 0.7])
xlabel(sp(2),'f (Hz)')
ylabel(sp(1),'S_x_x (m^2s^-^1)')
title(sp(1),'HTA, z/h_c = 0.03')
title(sp(2),'HTA, z/h_c = 0.33')
title(sp(3),'VTA, z/h_c = 0.04')
prettyfigures('text',13,'labels',14,'box',1)
% export_fig([figdir 'HTA_VTApsd_lsV2_1'],'-pdf')

