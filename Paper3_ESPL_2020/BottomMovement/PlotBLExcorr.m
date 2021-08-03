figure
sp(1) = subplot(211);
p(1) = plot(bdt,up);hold on
p(2) = plot(bdt,vp,'r');
p(3) = plot(bdt,wp,'g');
sp(2) = subplot(212);
plot(bdt,bdh,'k')
set(sp(1),'xticklabel',[])
set(sp,'xlim',[bdt(1) bdt(end)]);
datetickzoom('x','keepticks','keeplimits')
linkaxes(sp,'x')
leg = legend(p,{'u"';'v"';'w"'});
xlabel(sp(2),['Time on ' datestr(bdt(1),'dd-mm-yy')])
ylabel(sp(2),'BLE (m)')
ylabel(sp(1),'m/s'),title(sp(1),'Components of turbulence')
clear p

figure
sp(1) = subplot(211);
p(1) = plot(bdt,uw,'k');hold on
p(2) = plot(bdt,vw,'r');
sp(2) = subplot(212);
plot(bdt,bdh,'k')
set(sp(1),'xticklabel',[])
set(sp,'xlim',[bdt(1) bdt(end)]);
datetickzoom('x','keepticks','keeplimits')
linkaxes(sp,'x')
leg = legend(p,{'u"w"';'v"w"'});
xlabel(sp(2),['Time on ' datestr(bdt(1),'dd-mm-yy')])
ylabel(sp(2),'BLE (m)')
ylabel(sp(1),'m^2/s^2'),title(sp(1),'Reynolds stresses')

figure
sp(1) = subplot(211);
plot(bdt,k,'k')
sp(2) = subplot(212);
plot(bdt,bdh,'k')
set(sp(1),'xticklabel',[])
set(sp,'xlim',[bdt(1) bdt(end)]);
datetickzoom('x','keepticks','keeplimits')
linkaxes(sp,'x')
xlabel(sp(2),['Time on ' datestr(bdt(1),'dd-mm-yy')])
ylabel(sp(2),'BLE (m)')
ylabel(sp(1),'m^2/s^2'),title(sp(1),'TKE')

figure
sp(1) = subplot(211);
plot(bdt,ub,'k')
sp(2) = subplot(212);
plot(bdt,bdh,'k')
set(sp(1),'xticklabel',[])
set(sp,'xlim',[bdt(1) bdt(end)]);
datetickzoom('x','keepticks','keeplimits')
linkaxes(sp,'x')
xlabel(sp(2),['Time on ' datestr(bdt(1),'dd-mm-yy')])
ylabel(sp(2),'BLE (m)')
ylabel(sp(1),'m/s'),title(sp(1),'Near-bottom wave orbital velocity')

% export_fig(figure(4),[sdir expname '_UbBLE'],'-png')
% export_fig(figure(3),[sdir expname '_TKE_BLE'],'-png')
% export_fig(figure(2),[sdir expname '_RS_BLE'],'-png')
% export_fig(figure(1),[sdir expname '_turbBLE'],'-png')