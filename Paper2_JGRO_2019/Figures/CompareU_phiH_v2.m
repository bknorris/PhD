%Let's do some comparisons with other studies. For this script, I need to
%calculate Uc/Uinf from the Vectrino/Vector velocities (HTA) or lower/upper
%Vectrino velocities (VTA) and compare these to phi*(h/hc). Plot data on
%the same figure
%
%This is version 2 of this script
%
clear
load('d:\Mekong_W2015\DataAnalysis\Paper2\ADV_avs\CurrentHPhi_data.mat')
phiHhc = NaN(3,4);
phsd = NaN(3,4);
uwratio = NaN(3,4);
uwsd = NaN(3,4);
ucratio = NaN(3,4);
ucsd = NaN(3,4);
fn = fieldnames(mydata);
for i = 1:4
    dfn = fieldnames(mydata.(fn{i}));
    if i < 4
        for ii = 1:3
            phiHhc(ii,i) = nanmean(mydata.(fn{i}).(dfn{ii}).Hratio.*mydata.(fn{i}).vpro1.phi);
            phsd(ii,i) = nanstd(mydata.(fn{i}).(dfn{ii}).Hratio.*mydata.(fn{i}).vpro1.phi);
            uwratio(ii,i) = nanmean(mydata.(fn{i}).(dfn{ii}).Uwratio);
            uwsd(ii,i) = nanstd(mydata.(fn{i}).(dfn{ii}).Uwratio);
            ucratio(ii,i) = nanmean(mydata.(fn{i}).(dfn{ii}).Ucratio);
            ucsd(ii,i) = nanstd(mydata.(fn{i}).(dfn{ii}).Ucratio)/2;
        end
    else
        phiHhc(1,i) = nanmean(mydata.(fn{i}).Hratio.*mydata.(fn{i}).phi);
        phsd(1,i) = nanstd(mydata.(fn{i}).Hratio.*mydata.(fn{i}).phi);
        uwratio(1,i) = nanmean(mydata.(fn{i}).Uwratio);
        uwsd(1,i) = nanstd(mydata.(fn{i}).Uwratio);
        ucratio(1,i) = nanmean(mydata.(fn{i}).Ucratio);
        ucsd(1,i) = nanstd(mydata.(fn{i}).Ucratio)/2;
    end
end
uwratio(3,1) = 0.55;ucratio(3,1) = 0.11;
%Ok, now let's load other's data
fid = fopen('c:\Users\Bnorr\Documents\GradSchool\DataAnalysis\Paper2\StudyComparison\Compare.csv');
hdr = textscan(fgetl(fid),'%s','delimiter',',');hdr = hdr{:}';
cdata = textscan(fid,'%s%s%n%n%n','delimiter',',');
cdata = cell2struct(cdata,hdr,2);
figdir = 'c:\Users\Bnorr\Documents\GradSchool\DataAnalysis\Paper2\WorkingFigures\StudyComparison\';

%Plot the figure
f1 = figure(1);
set(f1,'PaperOrientation','portrait',...
    'position',[400 100 700 600]);
set(gcf,'color','w','paperpositionmode','auto')
cmap = csvread('c:\Users\Bnorr\Documents\GradSchool\Design\Qualitative8_3.csv');
cmap = cmap./255;
%Combine all the data for BF line
uwr = reshape(uwratio,12,1);uwr(isnan(uwr)) = [];
ucr = reshape(ucratio,12,1);ucr(isnan(ucr)) = [];
phr = reshape(phiHhc,12,1);phr(isnan(phr)) = [];
xs = [phr; phr; cdata.Phi.*cdata.hhc];
ys = [uwr; ucr; cdata.UcUinf];
%some nans in ys
nid = isnan(ys);ys(nid) = [];xs(nid) = [];
pv = polyfit(log(xs),log(ys),1);
y_hat=exp(pv(1)*log(xs)+pv(2));
%Do some statistics (R2, p-val)
yfit = polyval(pv,xs);
resid = ys-yfit;
SSresid = sum(resid.^2);
SStotal = (length(ys)-1)*var(ys);
rsq = 1-SSresid/SStotal;
[~,~,~,~,stats] = regress(ys,[ones(size(xs)) xs],0.05);
pval = stats(3);
disp(['r-squared: ' sprintf('%0.1f',rsq)])
disp(['p-value: ' sprintf('%0.4f',pval)])
bfl = loglog(xs,y_hat,'-k','linewidth',1.5);hold on
%Plot the other data
%Horstman2017 - low density
id = strcmp(cdata.CanopyType,'PLD');
xs = cdata.Phi(id).*cdata.hhc(id);
ys = cdata.UcUinf(id);
sp3 = plot(xs,ys,'^',...
    'color',cmap(1,:),...
    'linewidth',1.5,...
    'markerfacecolor','w',...
    'markersize',7);
%Horstman2017 - average density
id = strcmp(cdata.CanopyType,'PAD');
xs = cdata.Phi(id).*cdata.hhc(id);
ys = cdata.UcUinf(id);
sp4 = plot(xs,ys,'^',...
    'color',cmap(1,:),...
    'linewidth',1.5,...
    'markerfacecolor','w',...
    'markersize',7);
%Horstman2017 - high density
id = strcmp(cdata.CanopyType,'PHD');
xs = cdata.Phi(id).*cdata.hhc(id);
ys = cdata.UcUinf(id);
sp5 = plot(xs,ys,'^',...
    'color',cmap(1,:),...
    'linewidth',1.5,...
    'markerfacecolor','w',...
    'markersize',7);
%Horstman2017 - dowels, variable height
id = strcmp(cdata.CanopyType,'DowVar');
xs = cdata.Phi(id).*cdata.hhc(id);
ys = cdata.UcUinf(id);
sp6 = plot(xs,ys,'s',...
    'color',cmap(2,:),...
    'linewidth',1.5,...
    'markerfacecolor','w',...
    'markersize',7);
%Horstman2017 - dowels, uniform height
id = strcmp(cdata.CanopyType,'DowUni');
xs = cdata.Phi(id).*cdata.hhc(id);
ys = cdata.UcUinf(id);
sp7 = plot(xs,ys,'s',...
    'color',cmap(2,:),...
    'linewidth',1.5,...
    'markerfacecolor','w',...
    'markersize',7);
%Dunn 1996
id = strcmp(cdata.Study,'Dunn1996');
xs = cdata.Phi(id).*cdata.hhc(id);
ys = cdata.UcUinf(id);
sp8 = plot(xs,ys,'<',...
    'color',cmap(3,:),...
    'linewidth',1.5,...
    'markerfacecolor','w',...
    'markersize',7);
%Liu 2010
id = strcmp(cdata.Study,'Liu2010');
xs = cdata.Phi(id).*cdata.hhc(id);
ys = cdata.UcUinf(id);
sp9 = plot(xs,ys,'+',...
    'color',cmap(4,:),...
    'linewidth',1.5,...
    'markerfacecolor','w',...
    'markersize',7);
%G&N2004
id = strcmp(cdata.Study,'GN2004');
xs = cdata.Phi(id).*cdata.hhc(id);
ys = cdata.UcUinf(id);
sp10 = plot(xs,ys,'*',...
    'color',cmap(5,:),...
    'linewidth',1.5,...
    'markerfacecolor','w',...
    'markersize',7);
%Lowe2005
id = strcmp(cdata.CanopyType,'Unidirectional');
xs = cdata.Phi(id).*cdata.hhc(id);
ys = cdata.UcUinf(id);
sp11 = plot(xs,ys,'^',...
    'color',cmap(6,:),...
    'linewidth',1.5,...
    'markerfacecolor','w',...
    'markersize',7);
%Lowe2005
id = strcmp(cdata.CanopyType,'Wave_w');
xs = cdata.Phi(id).*cdata.hhc(id);
ys = cdata.UcUinf(id);
sp12 = plot(xs,ys,'o',...
    'color',cmap(6,:),...
    'linewidth',1.5,...
    'markerfacecolor',cmap(6,:),...
    'markersize',7);
%Lowe2005
id = strcmp(cdata.CanopyType,'Wave_c');
xs = cdata.Phi(id).*cdata.hhc(id);
ys = cdata.UcUinf(id);
sp13 = plot(xs,ys,'o',...
    'color',cmap(6,:),...
    'linewidth',1.5,...
    'markerfacecolor','w',...
    'markersize',7);
%Pujol2013
id = strcmp(cdata.CanopyType,'SRV_c');
xs = cdata.Phi(id).*cdata.hhc(id);
ys = cdata.UcUinf(id);
sp14 = plot(xs,ys,'p',...
    'color',cmap(7,:),...
    'linewidth',1.5,...
    'markerfacecolor','w',...
    'markersize',7);
%Pujol2013
id = strcmp(cdata.CanopyType,'SRV_w');
xs = cdata.Phi(id).*cdata.hhc(id);
ys = cdata.UcUinf(id);
sp15 = plot(xs,ys,'p',...
    'color',cmap(7,:),...
    'linewidth',1.5,...
    'markerfacecolor',cmap(7,:),...
    'markersize',7);
%Now plot my data
xs = phr;xerr = reshape(phsd,12,1);
xerr(isnan(xerr)) = [];
ys = uwr;yerr = reshape(uwsd,12,1);
yerr(isnan(yerr)) = [];
eb1 = ploterr(xs,ys,xerr,yerr,'.');
set(eb1,'linewidth',1.5,'color','k')
sp1 = plot(xs,ys,'d',...
    'color',cmap(8,:),...
    'linewidth',1.5,...
    'markerfacecolor',cmap(8,:),...
    'markersize',7);
%turn off the negative data warning
w = warning('query','last');
warning('off',w.identifier)
ys = ucr;yerr = reshape(ucsd,12,1);
yerr(isnan(yerr)) = [];
eb1 = ploterr(xs,ys,xerr,yerr,'.');
set(eb1,'linewidth',1.5,'color','k')
sp2 = plot(xs,ys,'d',...
    'color',cmap(8,:),...
    'linewidth',1.5,...
    'markerfacecolor','w',...
    'markersize',7);
%%%
prettyfigures('text',13,'labels',14,'box',1)
xlabel('\it\phi(h/h_c) [-]')
ylabel('\itU_c/U_\infty [-]')
set(gca,'xlim',[10^-4 10^0],...
    'ylim',[10^-3 10^1])
prettyfigures('text',13,'labels',14,'box',1)
bfleq = [sprintf('%0.2f',pv(1)) '(\it\phi(h/h_c))^{' sprintf('%0.2f',pv(2)) '}'];
legtxt = {'This Study';...
    'HM2017 - Pneum';'HM2017 - Dowels';'Dunn1996';'Liu2010';...
    'G-N2004';'LW2005 - Uni';'LW2005 - Osc';'PJ2013';bfleq};
leg = legend([sp2(1),sp4(1),sp6(1),sp8(1),...
    sp9(1),sp10(1),sp11(1),sp13(1),sp14(1),bfl],...
    legtxt,'fontsize',12,'linewidth',1.5);
set(leg,'position',[0.15 0.16 0.3 0.4])
xlabel('\it\phi(h/h_c) [-]')
ylabel('\itU_c/U_\infty [-]')
set(gca,'xlim',[10^-4 10^0],...
    'ylim',[0.01 2])
hold off
export_fig([figdir 'StudyComparison_v2'],'-pdf')
