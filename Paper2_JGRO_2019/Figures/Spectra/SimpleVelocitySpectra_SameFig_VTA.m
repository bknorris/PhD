%Simple velocity spectra version 2: this time, only plot x-shore
%velocities, all on top of each other on the same figure.
% clear
% load('C:\Users\bkn5\Projects\Mekong_W2015\DataAnalysis\Paper1\VTAvelocities.mat')
start = datenum(2015,03,14,07,50,00);stop = datenum(2015,03,14,07,52,00);
fn = fieldnames(VTA);
heading = 100;
inst = 3;
count = 1;
f1 = figure(1);
set(f1,'PaperOrientation','portrait',...
    'position',[400 100   800 800]);
set(gcf,'color','w','PaperPositionMode','auto')
c = [0 0 0;0.4 0.4 0.4;0.6 0.6 0.6];

for i = 2:4
    ind = find(VTA.(fn{i}).time >= start & VTA.(fn{i}).time <= stop);
    x = VTA.(fn{i}).x(ind,:);
    y = VTA.(fn{i}).y(ind,:);
    
    v = x(:,15); %just the middle bin
    u = y(:,15);
    
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
    nfft = floor(0.5*m);
    window = hanning(floor(0.05*m),'periodic');
    noverlap = length(0.7*length(window));
    [szu,fu] = pwelch(u,window,noverlap,nfft,fs);
    
    %windowing spectra: use 0.25Hz steps
    step = 0.25;
    nf = fs/2;
    ufits = zeros(length(fu),1);
    if i == 2
        window = 1;
        int = 1.5:step:8;
    else
        window  = 3;
        int = 1.5:step:10;
    end
    
    %calculate "moving best-fit line"
    for ii = 1:length(int)-window
        uix = find(fu >= int(ii) & fu <= int(ii+window));
        suss = mean(szu(uix));%fuss = fu(uix);
        %         pf = polyfit(log(fuss),log(suss),1);b = pf(1);yint = pf(2);
        %         y_hat = exp(b*log(fuss)+yint);
        %         ufits(uix) = y_hat;
        ufits(uix) = suss;
    end
    
    %fill low-f energy in fit
    uix = find(fu >= 0 & fu <= 1.6);
    ufits(uix) = szu(uix);
    
    %fill high-f energy in fit
    if i == 2
        uix = find(fu >= 8 & fu <= nf);
        ufits(uix) = szu(uix);
    else
        uix = find(fu >= 10 & fu <= nf);
        ufits(uix) = szu(uix);
    end
    
    %smooth for display
    uix = find(fu >= 1.5 & fu <= 18);
    ufits(uix) = smooth(ufits(uix),nf);
    
    P(count) = loglog(fu,ufits,'Color',c(count,:),'LineWidth',2); hold on
    count = count+1;
end
%-5/3 slope
int = 0.0022;
xs = linspace(2,45,length(szu));
ys = int.*(xs.^(-5/3));
loglog(xs,ys,'Color','r','LineWidth',1.5);
text(3.5,0.0005,'\itf ^-^5^/^3','FontName','Cambria',...
    'FontSize',14)
ylabel('S_u_u m^2s','FontName','Cambria',...
    'FontSize',14)
xlabel('Frequency (Hz)','FontName','Cambria',...
    'FontSize',14)
title('Cross-Shore ','FontName','Cambria',...
    'FontSize',14)
set(gca,'Xlim',[1E-2 1E2],'Ylim',[1E-6 1E0],...
    'LineWidth',1.5,'FontName','Cambria',...
    'FontSize',14)
leg = legend(P,{'z/h_c = 0.03';'z/h_c = 0.88';'z/h_c = 1.75'});set(leg,'Box','off')
hold off

export_fig(['c:\Users\bkn5\Projects\Mekong_W2015\Figures\Spectra&Waves\HorizSpectra\' 'VTAPowerSpectra'],'-pdf')