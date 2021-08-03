%Plot the elevation transect for better visualization of shore slope
%profile from mudflat -> forest. Data is from Sergio (e.g. Bryan et al.
%2017). 

clear
x = [-637.8183908,-635.2286677,-593.0727895,-553.7586377,-512.0008593,...
    -468.5035333,-426.695195,-384.6866411,-345.0467596,-307.844607,...
    -267.9753439,-226.3544854,-188.0661699,-150.0774648,-111.4211214,...
    -72.62019156,-37.34598617,-11.94479631,0];
z = [0.55565,0.529725,0.676125,0.535825,0.657825,0.76915,0.740175,0.75085,...
    0.75085,0.868275,0.8759,0.926225,0.868275,0.898775,0.968925,0.9857,...
    1.0772,1.17175,1.187];
xx = linspace(x(1),270,length(x));
zz = smooth(z,2);
% plot(xx,zz),hold on
x2 = xx(12:16);
z2 = zz(12:16);
% plot(x2,z2,'r')

f1 = figure(1);
set(f1,'PaperOrientation','portrait',...
    'position',[400 100   500 300]);
set(gcf,'color','w','PaperPositionMode','auto')
plot(xx,zz,'k','linewidth',1.5)
set(gca,...
    'xlim',[-80 100],...
    'xtick',-80:20:100,...
    'ylim',[0.8 1.2],...
    'ytick',0.8:0.2:1.2,...
    'linewidth',1.5,...
    'box','on',...
    'FontSize',14,...
    'FontName','Arial',...
    'TickDir','out')
xlabel('Cross Shore Distance (m)','fontsize',18)
ylabel('z (m)','fontsize',18)
figdir = 'd:\Projects\Mekong_W2015\Figures\Paper1\';
export_fig([figdir 'XShoreElevProf'],'-pdf','-nocrop')
