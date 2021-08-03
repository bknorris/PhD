fdir = 'f:\GradSchool\DataAnalysis\Paper2\WorkingFigures\ReynoldsStress\CF_V4\';
%End results: number of eliminated values
f1 = figure(1);
set(f1,'PaperOrientation','portrait',...
    'position',[400 100   600   400]);
set(gcf,'color','w','paperpositionmode','auto') 
imagesc(RS.day1.vpro3.time,0.06-0.04-linspace(0,0.03,35),RS.day1.vpro3.uw');
set(gca,'ydir','normal')
ylabel('Height above bed (m)')
datetick('x','HH:MM:SS','keepticks','keeplimits')
xlabel('Time on 07/03/15')
cb = colorbar;
ylabel(cb,'Reynolds Stress (Pa)')
title('VP3, x = 20 cm Reynolds Stress')
export_fig([fdir 'VP3_HTA1_RSestimates'],'-png')

%A good fit: Use vpro3, j = 46, jj = 4;
f2 = figure(2);
set(f2,'PaperOrientation','portrait',...
    'position',[400 100   800   600]);
set(gcf,'color','w','paperpositionmode','auto') 
subplot(211)
p1(1) = semilogx(k,oguwstar.*1.5,'k','linewidth',1.5);hold on
p1(2) = semilogx(k,oguw,'color',[0.6 0.6 0.6],'linewidth',1.5);
plot(ones(10,1)*kc,linspace(0,-6E-3,10),'color','k')
txt = '$\int{Co_{u''w''}dk} (m/s)$';
xlabel('k (rad/m)')
ylabel(txt,'interpreter','latex')
title('Ogive curves, a good fit: u_w/U < 2')

subplot(212)
p2(1) = semilogx(k,COuwvar,'color',[0.6 0.6 0.6],'linewidth',1.5);hold on
p2(2) = semilogx(k,COuwstar.*k*10,'k','linewidth',1.5);
plot(ones(10,1)*kc,linspace(-0.04,0.04,10),'color','k')
legend(p2,{'observed';'model'},'location','northeast')
txt = 'kCo_{u''w''} (m^2/s^2)';
xlabel('k (rad/m)')
ylabel(txt)
title('Variance-preserving cospectra, a good fit: u_w/U < 2')
% export_fig([fdir 'VP3_HTA1_goodfit'],'-png')

%A bad fit: use vpro3: 
f3 = figure(3);
set(f3,'PaperOrientation','portrait',...
    'position',[400 100   800   600]);
set(gcf,'color','w','paperpositionmode','auto') 
subplot(211)
p1(1) = semilogx(k,oguwstar,'k','linewidth',1.5);hold on
p1(2) = semilogx(k,oguw,'color',[0.6 0.6 0.6],'linewidth',1.5);
plot(ones(10,1)*kc,linspace(0,-0.12,10),'color','k')
set(gca,'ylim',[-0.12 0])
txt = '$\int{Co_{u''w''}dk} (m/s)$';
xlabel('k (rad/m)')
ylabel(txt,'interpreter','latex')
title('Ogive curves, a bad fit: u_w/U > 2')

subplot(212)
p2(1) = semilogx(k,COuwvar,'color',[0.6 0.6 0.6],'linewidth',1.5);hold on
p2(2) = semilogx(k,COuwstar.*k*5,'k','linewidth',1.5);
plot(ones(10,1)*kc,linspace(-0.3,0.04,10),'color','k')
set(gca,'ylim',[-0.3 0.05])
legend(p2,{'observed';'model'},'location','northeast')
txt = 'kCo_{u''w''} (m^2/s^2)';
xlabel('k (rad/m)')
ylabel(txt)
title('Variance-preserving cospectra, a bad fit: u_w/U > 2')
export_fig([fdir 'VP3_HTA1_poorfit'],'-png')
