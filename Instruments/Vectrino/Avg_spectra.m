%random plots: Spectral Averages
Xes = STAT.Spec.VelX.Spec+STAT.Spec.VelY.Spec;
Zes = (STAT.Spec.VelZ1.Spec+STAT.Spec.VelZ2.Spec)./2;
intv = 10;
f1 = figure;
set(f1,'PaperOrientation','portrait',...
    'position',[400 200   800   800]);
set(gcf,'color','w','PaperPositionMode','auto')
% suptitle('\bf\itS(X+Y) and S(Z1+Z2)/2')
[~,m] = size(Xes);
c = jet(m);
ax(1) = subplot(2,1,1);
for ii = 1:m
    area(STAT.Spec.f(:,ii),STAT.Spec.VelX.CI_up(:,ii),'FaceColor',[0.9 0.9 0.9],'LineStyle','none'),hold on
    area(STAT.Spec.f(:,ii),STAT.Spec.VelX.CI_low(:,ii),'FaceColor',[1 1 1],'LineStyle','none'),hold on
    p(ii) = plot(STAT.Spec.f(:,ii),Xes(:,ii),'-x','linewidth',1.5);
    set(p(ii),'Color',c(ii,:))
    hold on
    box on
%     xlabel('\bf\itHz')
    ylabel('\bf\itm^2/s')
%     title('\bf\itHorizontal Spectra S(X+Y)')
    names{ii} = sprintf('%d min to %d min',(intv*ii-intv),intv*ii);
    title('\bf\itS(X+Y)')
end
grid on
leg = legend(p,names,'location','northeast');
ax(2) = subplot(2,1,2);
for ii = 1:m
    area(STAT.Spec.f(:,ii),STAT.Spec.VelZ1.CI_up(:,ii),'FaceColor',[0.9 0.9 0.9],'LineStyle','none'),hold on
    area(STAT.Spec.f(:,ii),STAT.Spec.VelZ1.CI_low(:,ii),'FaceColor',[1 1 1],'LineStyle','none'),hold on
    p(ii) = plot(STAT.Spec.f(:,ii),Zes(:,ii),'-x','linewidth',1.5);
    set(p(ii),'Color',c(ii,:))
    hold on
    box on
    xlabel('\bf\itHz')
    ylabel('\bf\itm^2/s')
    title('\bf\itS(Z1+Z2)/2')
%     title('\bf\itAverage Vertical Spectra S(Z1+Z2)/2')
end
% leg = legend(p,names,'location','northeast');
grid on
set([ax(1) ax(2)],...
    'Xlimmode','auto',...
    'GridLineStyle',':')
set(ax(1),...
    'Xlim',[0 1],...
    'Ylim',[0 0.15],...
    'XTickLabel',[],...
    'YTick',0:0.05:0.15,...
    'position',[0.1 0.52 0.8 0.4])
set(ax(2),...
    'Xlim',[0 1],...
    'Ylim',[0 0.001],...
    'YTick',0:0.00025:0.001,...
    'position',[0.1 0.08 0.8 0.4])
prompt = 'Save Figure? [y/n] ';
result = input(prompt,'s');
if strcmp(result,'y') || strcmp(result,'yes');
    fpath = savefigdir;fname = [vname '_spectraAvgs'];
    export_fig([fpath fname],'-png','-m1','-r900','-opengl')
    disp(['Figure ' fname '.png saved'])
end
if strcmp(result,'n') || strcmp(result,'no');
end