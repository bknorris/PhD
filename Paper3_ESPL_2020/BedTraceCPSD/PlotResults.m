%Plot Results of h, vels, bd PSD and CPSD
%
% Paper 3: Sediment Motion in Mangroves
%
% This is Version 1.0 of this script
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

sdir = 'd:\Projects\Mekong_W2015\Figures\Paper3\BedCPSD\';

%Plot Spectra
sp = zeros(3,1);
f1 = figure(1);
set(f1,'PaperOrientation','portrait',...
    'position',[400 100   500  1000]);
sp(1) = subplot(311);
T = tb;
F = F2F2.(fn1{ii}).F;
Chh = F2F2.(fn1{ii}).Chh;
imagesc(T,F,10*log10(Chh)')
caxis([-100 -30])
cb = colorbar;
ylabel(cb,'dB')
datetickzoom('x','HH:MM','keepticks','keeplimits')

title('Depth')
sp(2) = subplot(312);
Cxx = F2F2.(fn1{ii}).Cxx;
imagesc(T,F,10*log10(Cxx)')
caxis([-100 -30])
cb = colorbar;
ylabel(cb,'dB')
datetickzoom('x','HH:MM','keepticks','keeplimits')

title('Cross-shore velocity')
sp(3) = subplot(313);
Cble = F2F2.(fn1{ii}).Cble;
imagesc(T,F,10*log10(Cble)')
caxis([-150 -100])
cb = colorbar;
ylabel(cb,'dB')
set(sp,'ylim',[0 1])
datetickzoom('x','HH:MM','keepticks','keeplimits')

xlabel(sp(3),['Time on ' datestr(vpt(1),'dd-mm-yy')])
ylabel(sp(1),'f (Hz)'),ylabel(sp(2),'f (Hz)')
ylabel(sp(3),'f (Hz)')
title('Bed Trace')
% prettyfigures('text',12,'labels',14,...
%     'box',1)

%Plot Example
start = datenum('11-Mar-2015 16:40:17');
T = F2F2.(fn1{ii}).time;
tmp = abs(T-start);
[~,idx] = min(tmp);
sp = zeros(3,1);
peaks = [0.0928 0.1714 0.2214 0.36];
pkid = zeros(4,1);
f2 = figure(2);
set(f2,'PaperOrientation','portrait',...
'position',[400 100   500  1000]);

sp(1) = subplot(311);
%depth
ff = F;
for i = 1:length(peaks)
tmp = abs(ff-peaks(i));
[~,ids] = min(tmp);
pkid(i) = ids;
end
F = F2F2.(fn1{ii}).F;
Chh = F2F2.(fn1{ii}).Chh(idx,:);
% L = mean(Chh-F2F2.aqdp2.CL(idx,:));
% U = mean(Chh+F2F2.aqdp2.CU(idx,:));
plot(F(pkid),Chh(pkid),'sq',...
'Color',[0.7 0.7 0.7],...
'MarkerSize',8,...
'MarkerFaceColor',[0.7 0.7 0.7]),hold on
plot(F,Chh,'-k',...
'LineWidth',1)
set(gca,'xscale','log','yscale','log')
% loglog(Fat,Chhc,'--r')
% errorbar(2,10E-3,L,U,'k')

sp(2) = subplot(312);
%velocities
Cxx = F2F2.(fn1{ii}).Cxx(idx,:);
Cyy = F2F2.(fn1{ii}).Cyy(idx,:);
plot(F(pkid),Cxx(pkid),'sq',...
'Color',[0.7 0.7 0.7],...
'MarkerSize',8,...
'MarkerFaceColor',[0.7 0.7 0.7]),hold on
plot(F,Cyy,'-',...
'Color',[0.5 0.5 0.5],...
'LineWidth',1)
plot(F,Cxx,'-k',...
'LineWidth',1)
set(gca,'xscale','log','yscale','log')

sp(3) = subplot(313);
%ble
Cble = F2F2.(fn1{ii}).Cble(idx,:);
plot(F(pkid),Cble(pkid),'sq',...
'Color',[0.7 0.7 0.7],...
'MarkerSize',8,...
'MarkerFaceColor',[0.7 0.7 0.7]),hold on
plot(F,Cble,'-k',...
'LineWidth',1)
set(gca,'xscale','log','yscale','log')
%global adjustments
set(sp,'xlim',[10^-2 1])
set(sp(3),'ylim',[10^-14 10^-8])
set([sp(1) sp(2)],'ylim',[10^-6 10^-2])
ylabel(sp(1),'Spectral Density (m^2Hz^-^1)')
ylabel(sp(2),'Spectral Density (m^2s^-^2Hz^-^1)')
ylabel(sp(3),'Spectral Density (m^2Hz^-^1)')
title(sp(1),'Depth'),title(sp(2),'Vectrino Velocity')
title(sp(3),'Bed Elevation')
xlabel(sp(1),'f (Hz)'),xlabel(sp(2),'f (Hz)'),xlabel(sp(3),'f (Hz)')

%
%
%
ci = cohere_signif_level(DOF);
f3 = figure(3);
set(f3,'PaperOrientation','portrait',...
    'position',[400 100   800 600]);
sp(1) = subplot(221);
MSChb = F2F2.(fn1{ii}).MSChb(idx,:);
plot(F(pkid),MSChb(pkid),'sq',...
    'Color',[0.7 0.7 0.7],...
    'MarkerSize',8,...
    'MarkerFaceColor',[0.7 0.7 0.7]),hold on
plot(F,ones(length(F),1)*ci,':',...
    'Color',[0.5 0.5 0.5])
plot(F,MSChb,'k')
set(gca,'xscale','log')

sp(2) = subplot(222);
MSCxb = F2F2.(fn1{ii}).MSCxb(idx,:);
plot(F(pkid),MSCxb(pkid),'sq',...
    'Color',[0.7 0.7 0.7],...
    'MarkerSize',8,...
    'MarkerFaceColor',[0.7 0.7 0.7]),hold on
plot(F,ones(length(F),1)*ci,':',...
    'Color',[0.5 0.5 0.5])
plot(F,MSCxb,'k')
set(gca,'xscale','log')

ciok = find(MSChb>ci);
sp(3) = subplot(223);
Phb = F2F2.(fn1{ii}).Phb(idx,:);
semilogx(F(ciok),Phb(ciok),'ok',...
    'linewidth',1.5)

ciok = find(MSCxb>ci);
sp(4) = subplot(224);
Pxb = F2F2.(fn1{ii}).Pxb(idx,:);
semilogx(F(ciok),Pxb(ciok),'ok',...
    'linewidth',1.5)
set(sp,'xlim',[10^-2 1])

%labeling
xlabel(sp(3),'f (Hz)')
xlabel(sp(4),'f (Hz)')
ylabel(sp(3),'Phase (degrees)')
ylabel(sp(1),'Coherence^2')
prettyfigures('text',12,'labels',14,...
    'box',1)
set(leg,'box','off')
