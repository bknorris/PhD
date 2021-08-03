%This script calculates the multiple linear regression coefficients of the
%turbulence/vegetation dataset. 

BigDataConglomerator

data = [reshape(eps,153,1) reshape(Hs,153,1) reshape(depth,153,1)...
    reshape(avg,153,1) reshape(zhc,153,1) reshape(X,153,1) reshape(phi,153,1)];

% Y = data(:,1);
% j=Y<0.005;
% X = [ones(length(data),1) data(:,2:end)];
% stepwise(X(j,:),Y(j))

%do some analysis with the data
%Stepwise reports that Epsilon is partially correlated with Hs, so
%normalize by Hs
Enorm = eps./Hs;En = reshape(Enorm,153,1);

%boxplot phi values on/off the mudflat
mp = data(:,7)==0;fp = data(:,7)~=0;
y1 = En(mp,:);y2 = En(fp,:);
bpx = NaN(length(En),2);
bpx(1:length(y1),1) = y1;
bpx(1:length(y2),2) = y2;bpx(bpx == 0) = NaN;

savefigdir = 'd:\Projects\Mekong_W2015\Figures\Paper2\';
f1 = figure(1);
set(f1,'PaperOrientation','portrait',...
    'position',[400 200   600   500]);
set(gcf,'color','w','PaperPositionMode','auto')
bp = boxplot(bpx,'colors','k','Labels',{'Mudflat','Fringe/Forest'});
set(bp,'linewidth',2);
set(gca,'Yscale','log','Ylim',[1E-6 1E-0],...
    'LineWidth',1.5,'FontSize',12,'FontName','Cambria')
ylabel('\epsilon*H_s^-^1 (ms^-^3)','FontSize',12,'FontName','Cambria')
set(findobj(gca,'Type','text'),'FontSize',12,'FontName','Cambria')
% export_fig([savefigdir 'EnviBoxPlot'],'-jpeg','-nocrop')

%also plot normalized epsilon by distance with phi and distance
Xs = [reshape(X,153,1);xqq];Phi = [reshape(phi,153,1); phiqq];
f2 = figure(2);
set(f2,'PaperOrientation','portrait',...
    'position',[400 200   1000   500]);
set(gcf,'color','w','PaperPositionMode','auto')
p(1) = plot(Xs,Phi,'o','MarkerEdgeColor','k',...
    'MarkerFaceColor','w','MarkerSize',10,...
    'LineWidth',1.5);hold on
%vertical line at zero
yl = linspace(0,1,100);xl = zeros(1,length(yl));
p(2) = plot(xl,yl,'--k','LineWidth',1.5);
axis([-80 100 0 5E-3])
ax1 = gca;
ax1_pos = get(ax1,'position'); % position of first axes
ax2 = axes('Position',ax1_pos,...
    'XAxisLocation','top',...
    'YAxisLocation','right',...
    'Color','none');
hold on
cmap = flipud([0,192,255; 0,129,255; 0,49,152]./256);
colormap(cmap)
for i = 1:3
    p(2) = errorbar(X(:,i),eps(:,i),stdev(:,i),'.',...
        'Color',cmap(i,:),'MarkerSize',8,'LineWidth',1.5,...
        'parent',ax2);hold on
end
set(ax2,'yscale','log','YLim',[1E-6 1E0],'xticklabel',[],'LineWidth',1.5,'box','on',...
    'FontSize',14,'FontName','Cambria')
ylabel(ax2,'\epsilon (m^2s^-^3)','FontSize',14,'FontName','Cambria')
xlabel(ax1,'Cross-Shore Distance (m)','FontSize',14,'FontName','Cambria')
ylabel(ax1,'\phi','FontSize',14,'FontName','Cambria')
set(ax1,'FontSize',14,'FontName','Cambria')
export_fig([savefigdir 'XshoreEandPhi'],'-jpeg','-nocrop')

%plot depth with distance for Karin
depth = reshape(depth,153,1);Xs = reshape(X,153,1);
f3 = figure(3);
set(f3,'PaperOrientation','portrait',...
    'position',[400 200   1000   500]);
set(gcf,'color','w','PaperPositionMode','auto')
p(1) = plot(Xs,depth,'+r','MarkerSize',10,...
    'LineWidth',1.5);hold on
%vertical line at zero
yl = linspace(0,2,100);xl = zeros(1,length(yl));
p(2) = plot(xl,yl,'--k','LineWidth',1.5);
axis([-80 100 0 1.6])
set(gca,'LineWidth',1.5,'box','on',...
    'FontSize',14,'FontName','Cambria')
xlabel('Cross-Shore Distance (m)','FontSize',14,'FontName','Cambria')
ylabel('Average Depth(m)','FontSize',14,'FontName','Cambria')
% export_fig([savefigdir 'XshoreDepth'],'-jpeg','-nocrop')

%plot wave height with distance for Karin
SigH = reshape(Hs,153,1);
f4 = figure(4);
set(f4,'PaperOrientation','portrait',...
    'position',[400 200   1000   500]);
set(gcf,'color','w','PaperPositionMode','auto')
p(1) = plot(Xs,SigH,'^r','MarkerSize',10,...
    'LineWidth',1.5);hold on
%vertical line at zero
yl = linspace(0,2,100);xl = zeros(1,length(yl));
p(2) = plot(xl,yl,'--k','LineWidth',1.5);
axis([-80 100 0 0.8])
set(gca,'LineWidth',1.5,'box','on',...
    'FontSize',14,'FontName','Cambria')
xlabel('Cross-Shore Distance (m)','FontSize',14,'FontName','Cambria')
ylabel('Significant Wave Height (m)','FontSize',14,'FontName','Cambria')
% export_fig([savefigdir 'XshoreHsig'],'-jpeg','-nocrop')

%plot wave height/depth with distance for Karin
HsD = SigH./depth;
f5 = figure(5);
set(f5,'PaperOrientation','portrait',...
    'position',[400 200   1000   500]);
set(gcf,'color','w','PaperPositionMode','auto')
p(1) = plot(Xs,HsD,'dr','MarkerSize',10,...
    'LineWidth',1.5);hold on
%vertical line at zero
yl = linspace(0,2,100);xl = zeros(1,length(yl));
p(2) = plot(xl,yl,'--k','LineWidth',1.5);
axis([-80 100 0 0.6])
set(gca,'LineWidth',1.5,'box','on',...
    'FontSize',14,'FontName','Cambria')
xlabel('Cross-Shore Distance (m)','FontSize',14,'FontName','Cambria')
ylabel('\gamma','FontSize',14,'FontName','Cambria')
% export_fig([savefigdir 'XshoreWavebreak'],'-jpeg','-nocrop')
