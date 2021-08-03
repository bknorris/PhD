clear, close all
fdir = 'g:\GradSchool\DataAnalysis\Paper2\WorkingFigures\ReynoldsStress\CF_V4\';
load('D:\Projects\Mekong_W2015\DataAnalysis\Paper2\RS_estimates\RStress_10min_v5_bdadj.mat')
%first, find all the variance preserving cospectra that are positive
%Cospectra
COuw = RS.day1.vpro3.COuw;
COuwstar = RS.day1.vpro3.COuwstar;
uw = RS.day1.vpro3.uw;
uwint = RS.day1.vpro3.uwint;
k = RS.day1.vpro3.k;
kCOuw = k.*COuw;
kCOuwstar = k.*COuwstar.*repmat(nanmean(uw,2),1,4097);
[m,n] = size(kCOuw);
pkCOuw = NaN(m,n);
pkCOuws = NaN(m,n);
pkoguw = NaN(m,n);
pkoguws = NaN(m,n);
%Ogive curves
oguw = RS.day1.vpro3.oguw;
oguwstar = cumsum(COuwstar.*repmat(nanmean(uw,2),1,4097).*0.4);
kk0 = NaN(m,n);
for j = 1:m
    %determine if kCOuw is negative or positive
    if abs(min(kCOuw(j,:))) < max(kCOuw(j,:))
        %maxima
        pkCOuw(j,:) = kCOuw(j,:)./nanmean(uw(j,:),2);
        pkCOuws(j,:) = kCOuwstar(j,:)./nanmean(uw(j,:),2);
        pkoguw(j,:) = oguw(j,:)./nanmean(uw(j,:),2);
        pkoguws(j,:) = oguwstar(j,:)./nanmean(uw(j,:),2);
        kk0(j,:) = k(j,:)./nanmean(RS.day1.vpro3.ku0(j,:),2);
    end
end
pkCOuw = nanmean(abs(pkCOuw));
pkCOuws = nanmean(pkCOuws);
pkoguw = nanmean(pkoguw);
pkoguws = nanmean(pkoguws);
kk0 = nanmean(kk0);

f1 = figure(1);
set(f1,'PaperOrientation','portrait',...
    'position',[400 100   800   400]);
set(gcf,'color','w','paperpositionmode','auto')
sp(1) = subplot(121);
bins = logspace(-2,2,30);
[b,~,q1,q3] = binmedian(kk0,(pkCOuw),bins);
h = ploterr(bins,b,[],{q1 q3},'o','logx','abshhx',0.01);hold on
set(h,...
    'markerfacecolor','k',...
    'color','k')
semilogx(kk0,pkCOuws,'-r','linewidth',1.5)
% set(gca,'xlim',[10^-2 10^2],'ylim',[-0.1 1.1])
xlabel('k/k_0')
ylabel('kCo_{uw}/{u''w''}')

sp(2) = subplot(122);
bins = logspace(-2,2,30);
[b,~,q1,q3] = binmedian(kk0,pkoguw,bins);
h = ploterr(bins,b,[],{q1 q3},'o','logx','abshhx',0.01);hold on
set(h,...
    'markerfacecolor','k',...
    'color','k')
semilogx(kk0,pkoguws,'-r','linewidth',1.5)
% set(gca,'xlim',[10^-2 10^2],'ylim',[-0.1 1.1],'yticklabel',[])
xlabel('k/k_0')
set(sp(1),'position',[0.15 0.15 0.35 0.75])
set(sp(2),'position',[0.57 0.15 0.35 0.75])
% export_fig([fdir 'VP3_HTA1_Cspc_ogiv'],'-pdf')

f2 = figure(2);
set(f2,'PaperOrientation','portrait',...
    'position',[400 100   500   400]);
set(f2,'color','w','paperpositionmode','auto')
x = nanmean(uwint,2);
y = nanmean(uw./1.5,2);
[x, y] = prepareCurveData(x,y);
pf = polyfit(x,y,1);
pv = polyval(pf,x);
rsq = rsquared(x,pv,length(pf));
plot(x,pv,'-k','linewidth',1.5),hold on
plot(x,y,'o','markersize',5,...
    'markerfacecolor',[0.6 0.6 0.6],...
    'color','k')
set(gca,'xlim',[-2E-3 3E-3],...
    'xtick',-2E-3:2E-3:4E-3,...
    'ylim',[-2E-3 3E-3],...
    'ytick',-2E-3:2E-3:4E-3,...
    'position',[0.15 0.15 0.8 0.78])
xlabel('Integral estimate of {u''w''} (m^2/s^2)')
ylabel('Fit estimate of {u''w''} (m^2/s^2)')
text(-1.8E-3,2E-3,sprintf('r^2 = %0.2f',rsq))
text(-1.8E-3,1.6E-3,sprintf('slope = %0.2f',pf(1)))
prettyfigures('text',12,'labels',14,'box',1)
export_fig([fdir 'VP3_HTA1_Int_Fit_uw_v2'],'-pdf')
