%plot AQDP and RBR pressure spectra from F2F2 for Karin.

rbr = load('RBRSpectrum.mat');
aqdp = load('AqdpSpectrum.mat');
sensors = catstruct(rbr.Spec,aqdp.Spec);
name = {'S80P';'HR3';'S81P';'AD5116';'S83P';'AD5117';'RD02'};
dfn = fieldnames(aqdp.Spec);
%need to plot the sensors in order from flats to forest - else use
%fieldnames
%PLOT%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
f1 = figure;
set(f1,'PaperOrientation','portrait',...
    'position',[400 200   1000   800]);
set(gcf,'color','w','PaperPositionMode','auto')
n = length(dfn);
c = linspecer(n,'qualitative');
ax = gca;
for ii = 1:n
    pq(ii) = plot(aqdp.Spec.(dfn{ii}).freq',aqdp.Spec.(dfn{ii}).psd,'-x','linewidth',1.5);
    set(pq(ii),'Color',c(ii,:))
    hold on
    box on  
end
leg = legend(pq,dfn,'northeast');
ylabel('\bf\itdbar^2/Hz','FontSize',14)
xlabel('\bf\itFrequency (Hz)','FontSize',14)
grid on
set(ax,'Xlim',[0 1],'Ylim',[0 25E-6],...
    'GridLineStyle',':',...
    'FontName', 'Helvetica','FontSize',14)

name = fieldnames(rbr.Spec);
%need to plot the sensors in order from flats to forest - else use
%fieldnames
%PLOT%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
f1 = figure;
set(f1,'PaperOrientation','portrait',...
    'position',[400 200   1000   800]);
set(gcf,'color','w','PaperPositionMode','auto')
n = length(name);
c = linspecer(n,'qualitative');
ax = gca;
for ii = 1:n
    pq(ii) = plot(rbr.Spec.(name{ii}).freq',rbr.Spec.(name{ii}).psd,'-x','linewidth',1.5);
    set(pq(ii),'Color',c(ii,:))
    hold on
    box on  
end
leg = legend(pq,name,'northeast');
ylabel('\bf\itdbar^2/Hz','FontSize',14)
xlabel('\bf\itFrequency (Hz)','FontSize',14)
grid on
set(ax,'Xlim',[0 1],'Ylim',[0 0.1],...
    'GridLineStyle',':',...
    'FontName', 'Helvetica','FontSize',14)