%Simple velocity spectra version 2: this time, only plot x-shore
%velocities, all on top of each other on the same figure.
clear
load('D:\Projects\Mekong_W2015\DataAnalysis\Paper1\HTAday1Vels.mat')
start = datenum(2015,03,07,14,24,00);stop = datenum(2015,03,07,14,26,00);
fn = fieldnames(HTA);
inst = 3;
count = 1;
f1 = figure();
set(f1,'PaperOrientation','portrait',...
    'position',[400 100   800 800]);
set(gcf,'color','w','PaperPositionMode','auto')
c = [0.6 0.6 0.6;0.4 0.4 0.4;0 0 0];

for i = 2:4
    ind = find(HTA.(fn{i}).time >= start & HTA.(fn{i}).time <= stop);
    x = HTA.(fn{i}).x(ind,:);
    y = HTA.(fn{i}).y(ind,:);
    
    v = x(:,9); %just the middle bin
    u = y(:,9);
    
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
    if i == 3
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
    if i == 3
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
leg = legend(P,{'x = -10cm';'x = 10cm';'x = 20cm'});set(leg,'Box','off')
hold off

% export_fig(['c:\Users\bkn5\Projects\Mekong_W2015\Figures\Spectra&Waves\HorizSpectra\' 'HTAday1PowerSpectra'],'-pdf')