%Let's do some comparisons with other studies. For this script, I need to
%calculate Uc/Uinf from the Vectrino/Vector velocities (HTA) or lower/upper
%Vectrino velocities (VTA) and compare these to phi*(h/hc)
%
%This is version 1 of this script
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
    'position',[400 100 1000 400]);
set(gcf,'color','w','paperpositionmode','auto')
sb(1) = subplot(121);
%Plot only the unidirectional data (i.e. all flume data) plus our
%currents/waves
id = strfind(cdata.CanopyType,'_w');
id = find(~not(cellfun('isempty',id)));
uwr = reshape(uwratio,12,1);uwr(isnan(uwr)) = [];
ucr = reshape(ucratio,12,1);ucr(isnan(ucr)) = [];
phr = reshape(phiHhc,12,1);phr(isnan(phr)) = [];
xs = [phr; phr; cdata.Phi(id).*cdata.hhc(id)];
ys = [uwr; ucr; cdata.UcUinf(id)];
%some nans in ys
nid = isnan(ys);ys(nid) = [];xs(nid) = [];
pv = polyfit(log(xs),log(ys),1);
y_hat=exp(pv(1)*log(xs)+pv(2));
bfl(1) = loglog(xs,y_hat,'-k','linewidth',1.5);hold on
%Now plot my data
xs = [phr; phr];xerr = [reshape(phsd,12,1); reshape(phsd,12,1)];
xerr(isnan(xerr)) = [];
ys = [uwr; ucr];yerr = [reshape(uwsd,12,1); reshape(ucsd,12,1)];
yerr(isnan(yerr)) = [];
eb1 = ploterr(xs,ys,xerr,yerr,'.k');
set(eb1,'linewidth',1.5)
sp1 = plot(xs,ys,'dk',...
    'linewidth',1.5,...
    'markerfacecolor','w',...
    'markersize',6);
%Plot the other data
%Horstman2017 - low density
id = strcmp(cdata.CanopyType,'PLD');
xs = cdata.Phi(id).*cdata.hhc(id);
ys = cdata.UcUinf(id);
sp2 = plot(xs,ys,'^r',...
    'linewidth',1.5,...
    'markerfacecolor','w',...
    'markersize',6);
%Horstman2017 - average density
id = strcmp(cdata.CanopyType,'PAD');
xs = cdata.Phi(id).*cdata.hhc(id);
ys = cdata.UcUinf(id);
sp3 = plot(xs,ys,'^b',...
    'linewidth',1.5,...
    'markerfacecolor','w',...
    'markersize',6);
%Horstman2017 - high density
id = strcmp(cdata.CanopyType,'PHD');
xs = cdata.Phi(id).*cdata.hhc(id);
ys = cdata.UcUinf(id);
sp4 = plot(xs,ys,'^g',...
    'linewidth',1.5,...
    'markerfacecolor','w',...
    'markersize',6);
%Horstman2017 - dowels, variable height
id = strcmp(cdata.CanopyType,'DowVar');
xs = cdata.Phi(id).*cdata.hhc(id);
ys = cdata.UcUinf(id);
sp5 = plot(xs,ys,'sm',...
    'linewidth',1.5,...
    'markerfacecolor','w',...
    'markersize',6);
%Horstman2017 - dowels, uniform height
id = strcmp(cdata.CanopyType,'DowUni');
xs = cdata.Phi(id).*cdata.hhc(id);
ys = cdata.UcUinf(id);
sp6 = plot(xs,ys,'sc',...
    'linewidth',1.5,...
    'markerfacecolor','w',...
    'markersize',6);
%Dunn 1996
id = strcmp(cdata.Study,'Dunn1996');
xs = cdata.Phi(id).*cdata.hhc(id);
ys = cdata.UcUinf(id);
sp7 = plot(xs,ys,'<k',...
    'linewidth',1.5,...
    'markerfacecolor','w',...
    'markersize',6);
%Liu 2010
id = strcmp(cdata.Study,'Liu2010');
xs = cdata.Phi(id).*cdata.hhc(id);
ys = cdata.UcUinf(id);
sp8 = plot(xs,ys,'pb',...
    'linewidth',1.5,...
    'markerfacecolor','w',...
    'markersize',6);
%G&N2004
id = strcmp(cdata.Study,'GN2004');
xs = cdata.Phi(id).*cdata.hhc(id);
ys = cdata.UcUinf(id);
sp9 = plot(xs,ys,'*r',...
    'linewidth',1.5,...
    'markerfacecolor','w',...
    'markersize',6);
%Lowe2005
id = strcmp(cdata.CanopyType,'Unidirectional');
xs = cdata.Phi(id).*cdata.hhc(id);
ys = cdata.UcUinf(id);
sp10 = plot(xs,ys,'oy',...
    'linewidth',1.5,...
    'markerfacecolor','y',...
    'markersize',6);

sb(2) = subplot(122);
%First, plot a best-fit line of wave data
%All wave data:
id = strfind(cdata.CanopyType,'_w');
id = find(not(cellfun('isempty',id)));
uwr = reshape(uwratio,12,1);uwr(isnan(uwr)) = [];
phr = reshape(phiHhc,12,1);phr(isnan(phr)) = [];
xs = [phr; cdata.Phi(id).*cdata.hhc(id)];
ys = [uwr; cdata.UcUinf(id)];
pv = polyfit(log(xs),log(ys),1);
y_hat=exp(pv(1)*log(xs)+pv(2));
bfl(2) = loglog(xs,y_hat,':k','linewidth',1.5);hold on
%Now plot my wave data
xs = phr;xerr = reshape(phsd,12,1);xerr(isnan(xerr)) = [];
ys = uwr;yerr = reshape(uwsd,12,1);yerr(isnan(yerr)) = [];
eb2 = ploterr(xs,ys,xerr,yerr,'.k');
set(eb2,'linewidth',1.5)
sp11 = plot(xs,ys,'dk',...
    'linewidth',1.5,...
    'markerfacecolor','w',...
    'markersize',6);
%Plot the other wave data
%Lowe2005
id = strcmp(cdata.CanopyType,'Wave_w');
xs = cdata.Phi(id).*cdata.hhc(id);
ys = cdata.UcUinf(id);
sp12 = plot(xs,ys,'^r',...
    'linewidth',1.5,...
    'markerfacecolor','w',...
    'markersize',6);
%Pujol2013
id = strcmp(cdata.CanopyType,'SRV_w');
xs = cdata.Phi(id).*cdata.hhc(id);
ys = cdata.UcUinf(id);
sp13 = plot(xs,ys,'og',...
    'linewidth',1.5,...
    'markerfacecolor','w',...
    'markersize',6);
%Now plot all current data:
id = strfind(cdata.CanopyType,'_c');
id = find(not(cellfun('isempty',id)));
ucr = reshape(ucratio,12,1);ucr(isnan(ucr)) = [];
xs = [phr; cdata.Phi(id).*cdata.hhc(id)];
ys = [ucr; cdata.UcUinf(id)];
pv = polyfit(log(xs),log(ys),1);
y_hat=exp(pv(1)*log(xs)+pv(2));
bfl(3) = loglog(xs,y_hat,'-k','linewidth',1.5);
%Now plot my wave data
xs = phr;xerr = reshape(phsd,12,1);xerr(isnan(xerr)) = [];
ys = ucr;yerr = reshape(ucsd,12,1);yerr(isnan(yerr)) = [];
eb3 = ploterr(xs,ys,xerr,yerr,'.k');
set(eb3,'linewidth',1.5)
sp14 = plot(xs,ys,'dk',...
    'linewidth',1.5,...
    'markerfacecolor','k',...
    'markersize',6);
%Plot the other current data
%Lowe2005
id = strcmp(cdata.CanopyType,'Wave_c');
xs = cdata.Phi(id).*cdata.hhc(id);
ys = cdata.UcUinf(id);
sp15 = plot(xs,ys,'^r',...
    'linewidth',1.5,...
    'markerfacecolor','r',...
    'markersize',6);
%Pujol2013
id = strcmp(cdata.CanopyType,'SRV_c');
xs = cdata.Phi(id).*cdata.hhc(id);
ys = cdata.UcUinf(id);
sp16 = plot(xs,ys,'og',...
    'linewidth',1.5,...
    'markerfacecolor','g',...
    'markersize',6);
prettyfigures('text',13,'labels',14,'box',1)
xlabel(sb(1),'\it\phi(h/h_c) [-]')
xlabel(sb(2),'\it\phi(h/h_c) [-]')
ylabel(sb(1),'\itUc/U_i_n_f [-]')
set(sb,'xlim',[10^-4 10^0],...
    'ylim',[10^-2 10^1])
set(sb(1),'position',[0.1 0.18 0.36 0.75])
set(sb(2),'position',[0.55 0.18 0.36 0.75],...
    'yticklabel',[])
% export_fig([figdir 'StudyComparison_v1'],'-pdf')
