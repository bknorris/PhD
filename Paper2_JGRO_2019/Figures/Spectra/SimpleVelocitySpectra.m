clear
load('D:Projects\Mekong_W2015\DataAnalysis\Paper1\HTAday2Vels.mat')

start = datenum(2015,03,08,15,20,00);stop = datenum(2015,03,08,15,26,00);
fn = fieldnames(HTA);
count = 1;
f1 = figure(count);
set(f1,'PaperOrientation','portrait',...
    'position',[400 100   1200 1000]);
set(gcf,'color','w','PaperPositionMode','auto')

for i = 2:4
ind = find(HTA.(fn{i}).time >= start & HTA.(fn{i}).time <= stop);
x = HTA.(fn{i}).x(ind,:);
y = HTA.(fn{i}).y(ind,:);

v = x(:,15); %just the middle bin
u = y(:,15);
% mu = nanmean(u);sigma = nanstd(u);
% pdfu = pdf('normal',u,mu,sigma);
% 
% %calc rms velocity
% urms = sqrt((sum(u.^2)/numel(u)));
% xu = u./urms;
% plot(xu,pdfu,'*')

%welch's method for psd
fs = 50;
[m,n] = size(u);
u = detrend(u);
v = detrend(v);
nfft = floor(0.5*m);
window = hanning(0.25*m,'periodic');
noverlap = 0.7*length(window);
[szu,fu] = pwelch(u,window,noverlap,nfft,fs);
[szv,fv] = pwelch(v,window,noverlap,nfft,fs);

%fit a line to various f intervals: check turbulent bandwidth?
%plot commands: figure 1
% figure(count)
% fmin = 3;fmax = 10;
% int = find(fu >= fmin & fu <= fmax);
% suss = szu(int);fuss = fu(int);
% pf = polyfit(log(fuss),log(suss),1);b = pf(1);yint = pf(2);
% y_hat=exp(b*log(fuss)+yint);
% %plot samples f(3,10Hz), linear fit for u velocity
% loglog(fuss,suss,'o','Color',[0.5 0.5 0.5]), hold on
% loglog(fuss,y_hat,'--r','LineWidth',1.5)
% fmin = 10;fmax = 25;
% int = find(fu >= fmin & fu <= fmax);
% suss = szu(int);fuss = fu(int);
% pf = polyfit(log(fuss),log(suss),1);b = pf(1);yint = pf(2);
% y_hat=exp(b*log(fuss)+yint);
% %plot samples f(10,25Hz), linear fit for u velocity
% loglog(fuss,suss,'o','Color',[0.5 0.5 0.5]), hold on
% loglog(fuss,y_hat,'--r','LineWidth',1.5)
int = 0.01;int2 = 0.1;
xs = linspace(3,25,length(szu));
xs2 = linspace(3,25,length(szu));
ys = int.*(xs.^(-5/3));
ys2 = int2.*(xs2.^(-5/2));
%plot -5/3 and -5/2 lines
% loglog(xs,ys,'Color','k','LineWidth',1.5);
% loglog(xs2,ys2,'--','Color','k','LineWidth',1.5);
% text(7,0.0005,'\itf^-^5^/^3')
% text(9,0.0005,'\itf^-^5^/^2')

%labeling commands: figure 1
% dt = datevec(stop-start);
% if dt(5) == 0
%     dt = dt(6); %# of seconds
%     title([fn{i} ': ' num2str(dt) ' second spectra, starting ' datestr(start,'dd/mm HH:MM')])
% else
%     dt = dt(5); %# of minutes
%     title([fn{i} ': ' num2str(dt) ' minute spectra, starting ' datestr(start,'dd/mm HH:MM')])
% end
% ylabel('PSD')
% xlabel('Frequency (Hz)')

%plot commands: figure 2
subplot(3,3,count)
g(1) = loglog(fu,szu,'-r'); hold on
g(2) = loglog(fv,szv,'--k');
loglog(xs,ys,'Color','k','LineWidth',1.5);
loglog(xs2,ys2,'--','Color','k','LineWidth',1.5);
text(10,0.001,'\itf^-^5^/^2')
text(30,0.0001,'\itf^-^5^/^3')
% set(gca,'Xlim',[10 25])
hold off
leg = legend([g(1),g(2)],{'S_u';'S_v'},'location','southwest');
set(gca,'Xlim',[1E-2 1E2],'Ylim',[1E-6 1E0])

%labeling commands: figure 2
title(fn{i})

ylabel('PSD')
xlabel('Frequency (Hz)')

%calculate linear fit of windows over short bandwidths to smooth the
%spectra. Look for the inflection point, it might tell you where to set the
%frequency window for the spectral slope figures!

%windowing spectra: use 1Hz steps with 5Hz window
window  = 3;
step = 1;
nf = fs/2;
ufits = zeros(length(fu),1);
vfits = zeros(length(fv),1);
int = 0:step:nf;

for ii = 1:length(int)-window
    uix = find(fu >= int(ii) & fu <= int(ii+window));
    suss = szu(uix);fuss = fu(uix);
    pf = polyfit(log(fuss),log(suss),1);b = pf(1);yint = pf(2);
    y_hat = exp(b*log(fuss)+yint);
    ufits(uix) = y_hat;

    vix = find(fv >= int(ii) & fv <= int(ii+window));
    svss = szv(vix);fvss = fv(vix);
    pf = polyfit(log(fvss),log(svss),1);b = pf(1);yint = pf(2);
    y_hat = exp(b*log(fvss)+yint);
    vfits(vix) = y_hat;
end
ufits = smooth(ufits,nf);
vfits = smooth(vfits,nf);

subplot(3,3,count+1)
loglog(fu,szu,'Color',[0.5 0.5 0.5]), hold on
loglog(fu,ufits,'-r','LineWidth',1.5)
loglog(xs,ys,'Color','k','LineWidth',1.5);
loglog(xs2,ys2,'--','Color','k','LineWidth',1.5);
text(10,0.001,'\itf^-^5^/^2')
text(30,0.0001,'\itf^-^5^/^3'), hold off
ylabel('PSD')
xlabel('Frequency (Hz)')
title([fn{i} ' Cross-Shore'])
set(gca,'Xlim',[1E-2 1E2],'Ylim',[1E-6 1E0])

subplot(3,3,count+2)
loglog(fv,szv,'Color',[0.5 0.5 0.5]), hold on
loglog(fv,vfits,'-r','LineWidth',1.5)
loglog(xs,ys,'Color','k','LineWidth',1.5);
loglog(xs2,ys2,'--','Color','k','LineWidth',1.5);
text(10,0.001,'\itf^-^5^/^2')
text(30,0.0001,'\itf^-^5^/^3'), hold off
ylabel('PSD')
xlabel('Frequency (Hz)')
title([fn{i} ' Along-Shore'])
set(gca,'Xlim',[1E-2 1E2],'Ylim',[1E-6 1E0])
count = count+3;
end
dt = datevec(stop-start);
if dt(5) == 0
    dt = dt(6); %# of seconds
    suptitle(['Spectral Analysis: ' num2str(dt) ' second spectra, starting ' datestr(start,'dd/mm HH:MM')])
else
    dt = dt(5); %# of minutes
    suptitle(['Spectral Analysis: ' num2str(dt) ' minute spectra, starting ' datestr(start,'dd/mm HH:MM')])
end
