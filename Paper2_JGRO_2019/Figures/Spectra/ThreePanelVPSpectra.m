%Simple velocity spectra version 3: this time, plot along, cross shore and
%vertical velocities, in three figures. We need to see how the shape of
%eddy structure changes in all three dimensions.

clear
load('C:\Users\bkn5\Projects\Mekong_W2015\DataAnalysis\Paper1\HTAday1Vels.mat')

start = datenum(2015,03,07,14,24,00);stop = datenum(2015,03,07,14,26,00);
fn = fieldnames(HTA);
inst = 3;
count = 1;
heading = 10;
f1 = figure(1);
set(f1,'PaperOrientation','landscape',...
    'position',[400 100 1000 500]);
set(gcf,'color','w','PaperPositionMode','auto')
c = [0.6 0.6 0.6;0.4 0.4 0.4;0 0 0];

for i = 2:4
    ind = find(HTA.(fn{i}).time >= start & HTA.(fn{i}).time <= stop);
    x = HTA.(fn{i}).x(ind,:);
    y = HTA.(fn{i}).y(ind,:);
    z = (HTA.(fn{i}).z1(ind,:)+HTA.(fn{i}).z2(ind,:))./2;
    v = x(:,2); %just the middle bin
    u = y(:,2);
    
    %rotate to the cross-shore direction
    rot = heading*pi/180;
    u = u.*(ones(size(u))*cos(rot)) + ...
        v.*(ones(size(v))*sin(rot));
    v = -u.*(ones(size(u))*sin(rot)) + ...
        v.*(ones(size(v))*cos(rot));
    
    %welch's method for psd
    fs = 50;
    [m,n] = size(u);
    u = detrend(u);
    v = detrend(v);
    z = detrend(z);
    nfft = floor(0.5*m);
    window = hanning(floor(0.05*m),'periodic');
    noverlap = length(0.7*length(window));
    [szu,fu] = pwelch(u,window,noverlap,nfft,fs);
    [szv,fv] = pwelch(v,window,noverlap,nfft,fs);
    [szz,fz] = pwelch(z,window,noverlap,nfft,fs);
    
    %windowing spectra: use 0.25Hz steps
    step = 0.25;
    nf = fs/2;
    ufits = zeros(length(fu),1);
    vfits = zeros(length(fv),1);
    zfits = zeros(length(fz),1);
    if i == 4
        window = 1;
        int = 1.5:step:8;
    else
        window  = 3;
        int = 1.5:step:10;
    end
   
    %calculate "moving best-fit line"
    for ii = 1:length(int)-window
        uix = find(fu >= int(ii) & fu <= int(ii+window));
        suss = mean(szu(uix));
        ufits(uix) = suss;
        
        vix = find(fv >= int(ii) & fv <= int(ii+window));
        svss = mean(szv(vix));
        vfits(vix) = svss;
        
        zix = find(fz >= int(ii) & fz <= int(ii+window));
        szss = mean(szz(zix));
        zfits(zix) = szss;
    end

    %fill low-f energy in fit
    uix = find(fu >= 0 & fu <= 1.6);
    ufits(uix) = szu(uix);
    
    vix = find(fv >= 0 & fv <= 1.6);
    vfits(vix) = szv(vix);
    
    zix = find(fz >= 0 & fz <= 1.6);
    zfits(zix) = szz(zix);
    
    %fill high-f energy in fit
    if i == 4
        uix = find(fu >= 8 & fu <= nf);
        ufits(uix) = szu(uix);
        
        vix = find(fv >= 8 & fv <= nf);
        vfits(vix) = szv(vix);
    
        zix = find(fz >= 8 & fz <= nf);
        zfits(zix) = szz(zix);
    else
        uix = find(fu >= 10 & fu <= nf);
        ufits(uix) = szu(uix);
        
        vix = find(fv >= 10 & fv <= nf);
        vfits(vix) = szv(vix);
    
        zix = find(fz >= 10 & fz <= nf);
        zfits(zix) = szz(zix);
    end
    
    %smooth for display
    uix = find(fu >= 1.5 & fu <= 11);
    ufits(uix) = smooth(ufits(uix),nf);
    
    vix = find(fu >= 1.5 & fu <= 11);
    vfits(vix) = smooth(vfits(vix),nf);
    
    zix = find(fu >= 1.5 & fu <= 11);
    zfits(zix) = smooth(zfits(zix),nf);
    
    %along shore
    pp(1) = subplot(131);
    P(count) = loglog(fv,vfits,'Color',c(count,:),'LineWidth',2); hold on
    %-5/3 slope
    int = 0.008;
    xs = linspace(3,50,length(szu));
    ys = int.*(xs.^(-5/3));
    loglog(xs,ys,'Color','r','LineWidth',1.5);
    text(6,0.0008,'\itf ^-^5^/^3','FontName','Cambria',...
        'FontSize',14)

    %cross shore
    pp(2) = subplot(132);
    R(count) = loglog(fu,ufits,'Color',c(count,:),'LineWidth',2); hold on
    %-5/3 slope
    loglog(xs,ys,'Color','r','LineWidth',1.5);
    text(6,0.0008,'\itf ^-^5^/^3','FontName','Cambria',...
        'FontSize',14)

    %vertical
    pp(3) = subplot(133);
    Q(count) = loglog(fz,zfits,'Color',c(count,:),'LineWidth',2); hold on
    %-5/3 slope
    int = 0.0022;
    xs = linspace(3,50,length(szu));
    ys = int.*(xs.^(-5/3));
    loglog(xs,ys,'Color','r','LineWidth',1.5);
    text(5,0.00045,'\itf ^-^5^/^3','FontName','Cambria',...
        'FontSize',14)
    count = count+1;
end

%position
set(pp(1),'position',[0.1 0.15 0.26 0.7])
set(pp(2),'position',[0.4 0.15 0.26 0.7],'YTickLabel',[])
set(pp(3),'position',[0.7 0.15 0.26 0.7],'YTickLabel',[])

%labeling
set(pp,'LineWidth',1.5,'FontName','Cambria',...
    'FontSize',14)
xlabel(pp(2),'Frequency (Hz)','FontName','Cambria',...
    'FontSize',14)
ylabel(pp(1),'PSD m^2s','FontName','Cambria',...
    'FontSize',14)

leg = legend(P,{'x = -10cm';'x = 10cm';'x = 20cm'});
set(leg,'position',[0.76 0.25 0.05 0.05],...
'Box','off','LineWidth',1.5)
hold off
export_fig(['c:\Users\bkn5\Projects\Mekong_W2015\Figures\Spectra&Waves\HorizSpectra\' 'HTAday1PSpec_xyz'],'-jpeg','-nocrop')