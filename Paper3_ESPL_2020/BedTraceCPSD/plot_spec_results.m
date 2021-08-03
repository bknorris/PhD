function plot_spec_results(expname,vars)
%Plot Results of h, vels, bd PSD and CPSD
%
% Paper 3: Sediment Motion in Mangroves
%
% This is Version 1.0 of this script
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

sdir = 'd:\Projects\Mekong_W2015\Figures\Paper3\BedCPSD\';
close all

%Plot Spectra
sp = zeros(3,1);
f1 = figure(1);
set(f1,'PaperOrientation','portrait',...
    'position',[400 100   500  1000]);
%DEPTH
sp(1) = subplot(311);
T = vars.h.time;
F = vars.h.f;
psd = vars.h.psd;
imagesc(T,F,10*log10(psd)')
caxis([-100 0])
cb = colorbar;
ylabel(cb,'dB')
datetickzoom('x','HH:MM','keepticks','keeplimits')
title('Depth')

%XSHORE VELOCITY
sp(2) = subplot(312);
psd = vars.x.psd;
imagesc(T,F,10*log10(psd)')
caxis([-100 0])
cb = colorbar;
ylabel(cb,'dB')
datetickzoom('x','HH:MM','keepticks','keeplimits')
title('Cross-shore velocity')

%BOTTOM TRACE
sp(3) = subplot(313);
psd = vars.bd.psd;
imagesc(T,F,10*log10(psd)')
caxis([-150 -80])
cb = colorbar;
ylabel(cb,'dB')
set(sp,'ylim',[0 1])
datetickzoom('x','HH:MM','keepticks','keeplimits')
title('Bed Trace')
xlabel(sp(3),['Time on ' datestr(T(1),'dd-mm-yy')])
ylabel(sp(1),'f (Hz)'),ylabel(sp(2),'f (Hz)')
ylabel(sp(3),'f (Hz)')
% prettyfigures('text',12,'labels',14,...
%     'box',1)

%Pick a point on Figure 1 to plot other figs from
disp('USER pick (1) point in t-s to plot other figures from')
disp('Use left click to zoom, right click to mark points')
disp('press RETURN when finished')
[xp,~] = ginput(1);


%Plot Example - PSD
start = datenum(xp);
T = vars.h.time;
tmp = abs(T-start);
[~,idx] = min(tmp);
%Use ginput to select peaks in spectra
yes = 0;
while yes == 0
    %depth
    f2 = figure(2);
    F = vars.h.f;
    psd = vars.h.psd(idx,:);
    
    ci = log10(vars.h.ci(idx,:));cx = ci(2)-ci(1);
    plot(F,psd,'-k',...
        'LineWidth',1)
    set(gca,'xscale','log','yscale','log'),hold on
    plot(0.8,0.02,'.k'), hold on
    sc = 5^-2;
    plot([0.8 0.8],[0.02-(cx/2)*sc 0.02+(cx/2)*sc],'-k')
    set(gca,'xlim',[10^-2 1])
    disp('USER pick points for spectral peaks')
    disp('Use left click to zoom, right click to mark points')
    disp('press RETURN when finished')
    
    peaks1 = zeros(length(xp),2);
    [xp,yp] = ginput_zoom;
    for i = 1:length(xp)
        tmp = abs(F-xp(i));
        [~,idxx] = min(tmp);
        tmp = abs(psd-yp(i));
        [~,idyy] = min(tmp);
        peaks1(i,1) = idxx;
        peaks1(i,2) = idyy;
    end
    plot(F(peaks1(:,1)),psd(peaks1(:,2)),'sq',...
        'color',[0.7 0.7 0.7],...
        'markerfacecolor',[0.7 0.7 0.7],...
        'markersize',8)
    
    f3 = figure(3);
    F = vars.x.f;
    psd = vars.x.psd(idx,:);
    
    ci = log10(vars.h.ci(idx,:));cx = ci(2)-ci(1);
    plot(F,psd,'-k',...
        'LineWidth',1)
    set(gca,'xscale','log','yscale','log'),hold on
    plot(0.8,0.02,'.k')
    sc = 5^-2;
    plot([0.8 0.8],[0.02-(cx/2)*sc 0.02+(cx/2)*sc],'-k')
    set(gca,'xlim',[10^-2 1])

    disp('USER pick points for spectral peaks')
    disp('Use left click to zoom, right click to mark points')
    disp('press RETURN when finished')
    
    
    [xp,yp] = ginput_zoom;
    peaks2 = zeros(length(xp),2);
    for i = 1:length(xp)
        tmp = abs(F-xp(i));
        [~,idxx] = min(tmp);
        tmp = abs(psd-yp(i));
        [~,idyy] = min(tmp);
        peaks2(i,1) = idxx;
        peaks2(i,2) = idyy;
    end
    plot(F(peaks2(:,1)),psd(peaks2(:,2)),'d',...
        'color',[0.7 0.7 0.7],...
        'markerfacecolor',[0.7 0.7 0.7],...
        'markersize',8)
    
    f4 = figure(4);
    F = vars.bd.f;
    psd = vars.bd.psd(idx,:);
    ci = log10(vars.bd.ci(idx,:));cx = ci(2)-ci(1);
    plot(F(peaks1(:,1)),psd(peaks1(:,1)),'sq',...
        'Color',[0.7 0.7 0.7],...
        'MarkerSize',8,...
        'MarkerFaceColor',[0.7 0.7 0.7]),hold on
    plot(F(peaks2(:,1)),psd(peaks2(:,1)),'d',...
        'Color',[0.7 0.7 0.7],...
        'MarkerSize',8,...
        'MarkerFaceColor',[0.7 0.7 0.7])
    plot(F,psd,'-k',...
        'LineWidth',1)
    set(gca,'xscale','log','yscale','log')
    sc = 5^-14;
    plot(0.6,10^-10,'.k')
    plot([0.6 0.6],[10^-10-(cx/2)*sc 10^-10+(cx/2)*sc],'-k')
    set(gca,'xlim',[10^-2 1])
    
    prompt = 'Is the selection correct [y/n]? ';
    result = input(prompt,'s');

    if strcmp(result,'y');
        close(f2,f3,f4),yes = 1;
    elseif strcmp(result,'n');
        close(f2,f3,f4)
        continue
    end
end

f2 = figure(2);
set(f2,'PaperOrientation','portrait',...
    'position',[400 100   500  1000]);

sp(1) = subplot(311);
%depth
F = vars.h.f;
psd = vars.h.psd(idx,:);
ci = log10(vars.h.ci(idx,:));cx = ci(2)-ci(1);
plot(F(peaks1(:,1)),psd(peaks1(:,1)),'sq',...
    'Color',[0.7 0.7 0.7],...
    'MarkerSize',8,...
    'MarkerFaceColor',[0.7 0.7 0.7]),hold on
plot(F(peaks2(:,1)),psd(peaks2(:,1)),'d',...
    'Color',[0.7 0.7 0.7],...
    'MarkerSize',8,...
    'MarkerFaceColor',[0.7 0.7 0.7])
plot(F,psd,'-k',...
    'LineWidth',1)
set(gca,'xscale','log','yscale','log')
sc = 5^-2;
plot(0.6,0.02,'.k')
plot([0.6 0.6],[0.02-(cx/2)*sc 0.02+(cx/2)*sc],'-k')
text(0.65,0.02,'95%')
text(2E-2,5E-4,'a) Depth')
set(gca,'xlim',[10^-2 1])

sp(2) = subplot(312);
%depth
F = vars.x.f;
psd1 = vars.x.psd(idx,:);
psd2 = vars.y.psd(idx,:);
ci = log10(vars.h.ci(idx,:));cx = ci(2)-ci(1);
plot(F(peaks1(:,1)),psd1(peaks1(:,1)),'sq',...
    'Color',[0.7 0.7 0.7],...
    'MarkerSize',8,...
    'MarkerFaceColor',[0.7 0.7 0.7]),hold on
plot(F(peaks2(:,1)),psd1(peaks2(:,1)),'d',...
    'Color',[0.7 0.7 0.7],...
    'MarkerSize',8,...
    'MarkerFaceColor',[0.7 0.7 0.7])
lg = zeros(2,1);
lg(1) = plot(F,psd1,'-k',...
    'LineWidth',1);
lg(2) = plot(F,psd2,'color',[0.5 0.5 0.5],...
    'LineWidth',1);
set(gca,'xscale','log','yscale','log')
sc = 0.5;
plot(0.6,0.2,'.k')
plot([0.6 0.6],[0.2-(cx/2)*sc 0.2+(cx/2)*sc],'-k')
text(0.65,0.2,'95%')
text(2E-2,8E-4,'b) Velocity')
leg = legend(lg,{'Cross-shore';'Along-shore'});
set(leg,'position',[0.6 0.42 0.05 0.05])
set(gca,'xlim',[10^-2 1])%

sp(3) = subplot(313);
%depth
F = vars.bd.f;
psd = vars.bd.psd(idx,:);
ci = log10(vars.bd.ci(idx,:));cx = ci(2)-ci(1);
plot(F(peaks1(:,1)),psd(peaks1(:,1)),'sq',...
    'Color',[0.7 0.7 0.7],...
    'MarkerSize',8,...
    'MarkerFaceColor',[0.7 0.7 0.7]),hold on
plot(F(peaks2(:,1)),psd(peaks2(:,1)),'d',...
    'Color',[0.7 0.7 0.7],...
    'MarkerSize',8,...
    'MarkerFaceColor',[0.7 0.7 0.7])
plot(F,psd,'-k',...
    'LineWidth',1)
set(gca,'xscale','log','yscale','log')
sc = 5^-14;
plot(0.6,10^-10,'.k')
plot([0.6 0.6],[10^-10-(cx/2)*sc 10^-10+(cx/2)*sc],'-k')
text(0.65,10^-10,'95%')
text(2E-2,1E-11,'c) BLE')
set(gca,'xlim',[10^-2 1])

%global adjustments
set(sp,'xlim',[10^-2 1])
set(sp(3),'ylim',[10^-12 10^-8])
set(sp(1),'ylim',[10^-4 10^-1])
set(sp(2),'ylim',[10^-4 10^-0])
ylabel(sp(1),'Spectral Density (m^2Hz^-^1)')
ylabel(sp(2),'Spectral Density (m^2s^-^2Hz^-^1)')
ylabel(sp(3),'Spectral Density (m^2Hz^-^1)')
xlabel(sp(1),'f (Hz)'),xlabel(sp(2),'f (Hz)'),xlabel(sp(3),'f (Hz)')
% 
f3 = figure(3);
set(f3,'PaperOrientation','portrait',...
    'position',[400 100   800 600]);
sp(1) = subplot(221);
F = vars.hbd.f;
msc = vars.hbd.coh(idx,:);
ci = vars.hbd.ci(idx,:);
plot(F(peaks1(:,1)),msc(peaks1(:,1)),'sq',...
    'Color',[0.7 0.7 0.7],...
    'MarkerSize',8,...
    'MarkerFaceColor',[0.7 0.7 0.7]),hold on
plot(F(peaks2(:,1)),msc(peaks2(:,1)),'d',...
    'Color',[0.7 0.7 0.7],...
    'MarkerSize',8,...
    'MarkerFaceColor',[0.7 0.7 0.7]),hold on
plot(F,ones(length(F),1)*ci,':',...
    'Color','k')
plot(F,msc,'k')
set(gca,'xscale','log')

ciok = find(msc>ci);
sp(3) = subplot(223);
phase = vars.hbd.phase(idx,:);
semilogx(F(ciok),phase(ciok),'ok',...
    'linewidth',1.5)

sp(2) = subplot(222);
msc = vars.xbd.coh(idx,:);
ci = vars.xbd.ci(idx,:);
plot(F(peaks1(:,1)),msc(peaks1(:,1)),'sq',...
    'Color',[0.7 0.7 0.7],...
    'MarkerSize',8,...
    'MarkerFaceColor',[0.7 0.7 0.7]),hold on
plot(F(peaks2(:,1)),msc(peaks2(:,1)),'d',...
    'Color',[0.7 0.7 0.7],...
    'MarkerSize',8,...
    'MarkerFaceColor',[0.7 0.7 0.7]),hold on
plot(F,ones(length(F),1)*ci,':',...
    'Color','k')
plot(F,msc,'k')
set(gca,'xscale','log')

ciok = find(msc>ci);
sp(4) = subplot(224);
phase = vars.hbd.phase(idx,:);
semilogx(F(ciok),phase(ciok),'ok',...
    'linewidth',1.5)
set(sp,'xlim',[10^-2 1])

%labeling
title(sp(1),'Depth-BLE')
title(sp(2),'Velocity-BLE')
xlabel(sp(3),'f (Hz)')
xlabel(sp(4),'f (Hz)')
ylabel(sp(3),'Phase (degrees)')
ylabel(sp(1),'Coherence^2')
prettyfigures('text',12,'labels',12,...
    'box',1)
set(leg,'box','off')

fname = expname;
% export_fig(f1,[sdir fname 'hVbd_PSD_color'],'-png')
% export_fig(f2,[sdir fname 'hVbd_PSD_exmp'],'-png')
% export_fig(f3,[sdir fname 'hVbd_CPSD_exmp'],'-png')
end