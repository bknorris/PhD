clear, close all
f1 = figure(1);
set(f1,'PaperOrientation','portrait',...
    'position',[400 100   500   400],...
    'renderer','painters');hold on

load('e:\Mekong_W2015\DataAnalysis\Paper3\BedStress\06-03-15\F2F2_06_vpro1.mat')
xs = dat.xfit(9064,:);ys = dat.yfit(9064,:);
id = find(isnan(xs),1,'first');
xs(id:end) = linspace(xs(id-1),0,length(id:35));
ys(id:end) = linspace(ys(id-1),-7,length(id:35));
%line of best fit
xi = xs(1:id)./xs(id);yi = ys(1:id);
pf = polyfit(xi,yi,1);
xx = linspace(0.74,1,10);
pv = polyval(pf,xx);
%highlight points of BFL
plot(xi,yi,'.',...
    'markersize',7,...
    'color','k'),hold on
plot(xx,pv,'-r')
set(gca,...
    'xlim',[0.74 1.02],...
    'ylim',[-5.2, -3.4])

title('March 6^{th}, 2015')
ylabel('ln(z) [m]')
xlabel('u_z/u_{z = 0.02} [-]')
prettyfigures('text',12,'labels',13,'box',1,'grid',1,'gstyle',':')
sfdir = 'c:\Users\Bnorr\Documents\GradSchool\Writing\GM_2018_Norris\Figures\V6\';
% export_fig([sfdir 'z0_fit_v2'],'-pdf','-nocrop')
