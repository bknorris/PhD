%Plot MSC and velocity spectra to illustrate Bricker & Monismiths 2007's method
%for determining the wave band in spectra.

clear, close all
load('d:\Projects\Mekong_W2015\DataAnalysis\Paper2\Spectra\CutoffFreqSpectra.mat')
f1 = figure(1);
set(f1,'PaperOrientation','portrait',...
    'position',[400 100   400   500]);
set(gcf,'color','w','paperpositionmode','auto')   
sup(1) = subplot(211);
chl = coheresig;
semilogx(f,surfp_msc,'k',...
    'linewidth',1.5), hold on
plot(f,repmat(chl,length(f),1),'--r',...
    'linewidth',1.5)
f2 = f(surfp_msc>chl);
sp = surfp_msc(surfp_msc>chl);
fid = find(f2>=0.05&f2<=2);
x = f2(fid);y = sp(fid);
markx = x(1:20:end);marky = y(1:20:end);
plot(markx,marky,'o',...
    'linewidth',1.5,...
    'color',[0.6 0.6 0.6],...
    'markersize',6)
%
sup(2) = subplot(212);
loglog(f,Cww,'k',...
    'linewidth',1.5), hold on
cw = Cww(surfp_msc>chl);
y = cw(fid);
marky = y(1:20:end);
plot(markx,marky,'o',...
    'linewidth',1.5,...
    'color',[0.6 0.6 0.6],...
    'markersize',6)
%turbulent band line fit
fc = find(f >= 1.4 & f <= 8);
inert = Cww(fc);ff = f(fc);
p = polyfit(log(ff),log(inert),1);
y_hat=exp(p(1)*log(f)+p(2));
loglog(f(fc),y_hat(fc),'-r','linewidth',2)
%-5/3 slope line
xs = linspace(2,40,length(Cww));
ys = 2E-3.*(xs.^(-5/3));
plot(xs,ys,'-.k','LineWidth',1.5);hold on
text(10,2E-4,'f^{ -5/3}')
%Global plot adjustments
set(sup,'xlim',[10^-2 10^2])
set(sup(1),'ylim',[0 1],...
    'ytick',0:0.5:1,...
    'xticklabel',[],...
    'position',[0.19 0.62 0.77 0.3])
set(sup(2),'ylim',[10^-6 10^0],...
    'position',[0.19 0.15 0.77 0.4])
xlabel(sup(2),'f (Hz)')
ylabel(sup(2),'S_{ww} (m^2/s^2/Hz)')
ylabel(sup(1),'Coherence Squared')
prettyfigures('text',12,'labels',14,'box',1)
figdir = 'f:\GradSchool\DataAnalysis\Paper2\WorkingFigures\Spectra\';
export_fig([figdir 'BM_2007_CutoffFreq'],'-pdf')

