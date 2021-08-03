%plot Phi from cyl and phi from vol
% p1 = 'D:\Projects\Code\ImageProcessing\';
% run([p1 'CatVegGeometry.m'])
vpdir = 'd:\Projects\Documents\Writing\DataReports\';
files = dir([vpdir '*_vol.mat']);files = {files.name};

% phic = [vegdat.four.Phic(1:5,:); vegdat.four.Phic(16:end,:); vegdat.five.Phic];
% phi = [vegdat.four.Phi(1:5,:); vegdat.four.Phi(16:end,:); vegdat.five.Phi];
phic = zeros(54,4);phi = zeros(54,4);
for i = 1:4
    load([vpdir files{i}])
    phic(:,i) = vegdat.phi_cyl;
    phi(:,i) = vegdat.phi_vol;
end
phi = mean(phi,2);phic = mean(phic,2);
%linear fit
pf = polyfit(phic,phi,1);
xs = linspace(0,1,length(phic));
pv = polyval(pf,xs);

%figure
f1 = figure(1);
set(f1,'PaperOrientation','portrait',...
    'position',[400 200   600   500]);
set(gcf,'color','w','PaperPositionMode','auto')
plot(xs,pv,'LineWidth',1.5,...
    'Color','k')
hold on
plot(phic,phi,'o',...
    'LineWidth',1.5,...
    'MarkerSize',10,...
    'MarkerFaceColor',[0.8 0.8 0.8],...
    'MarkerEdgeColor','k')

axis([0 0.7 0 0.04])
set(gca,'Box','on',...
    'LineWidth',1.5,...
    'TickDir','out',...
    'FontName','Arial',...
    'FontSize',13)
xlabel('\phi_c_y_l','FontSize',15)
ylabel('\phi_v_o_l','FontSize',15)
%display slope on plot
text(0.2,0.032,['\phi_v_o_l / \phi_c_y_l = ' sprintf('%0.2f',pf(1))],'FontSize',14)
figdir = 'd:\Projects\Mekong_W2015\Figures\Paper2\';
export_fig([figdir 'Phi_cyl_phi_vol_20cm_avg'],'-jpeg','-nocrop')