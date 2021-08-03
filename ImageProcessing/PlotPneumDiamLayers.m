%plot quadrat layers

% load('C:\users\bkn5\Projects\Mekong_F2014\Images\Quadrats\Reconstructions\Q4\Analysis\NoQuad\Q4stats_fixed.mat')
ddir = 'C:/Users/bkn5/Projects/Mekong_W2015/Images/Quadrats/Q2/Reconstructions/2B/Analysis/';
fdir = 'C:/Users/bkn5/Projects/Mekong_W2015/Figures/Quadrats/';
load([ddir 'Q2Bstats_fixed.mat'])

h = 0.005; %spacing between layers (m)
fn = fieldnames(DATA);n = length(fn);
H = 0:h:h*n; %vector of heights to max H
maxh = str2double(sprintf('%.2f',H(end)+0.1)); %for setting the upper z limit of the figure
%calculate categories of diameter based on the max(Radii) of bottom
%(largest) and top (smallest) layers in the model
maxd = mean(DATA.(fn{2}).Rradii(1,:)); %fn{1} is the info field of DATA
mind = mean(DATA.(fn{end}).Rradii(end,:));
diams = linspace(mind,maxd,n);

%(x,y,z) points of VecPros in the Quadrat, measured in CloudCompare
%[VP1 VP2 VP3]
vpx = [0.8033 0.5997 0.4978;0.8033 0.5997 0.4978;0.8033 0.5997 0.4978];
vpy = [0.245 0.245 0.245;0.245 0.245 0.245;0.254 0.254 0.254;];
vpz = [0.062 0.063 0.061;0.25 0.25 0.25;0.46 0.46 0.46]; %relative to the top of the quadrat, where the model bed level is initialized

f1 = figure;
set(f1,'PaperOrientation','portrait',...
    'position',[400 200   1000   800]);
set(gcf,'color','w','PaperPositionMode','auto')
c = winter(n);colormap(winter)
%generate 'ground' and borders in figure
patch([0 1 1 0],[0 0 1 1],[0.85 0.85 0.85])
patch([0 0 0 0],[1 0 0 1],[0 0 maxh maxh],[0.85 0.85 0.85])
patch([1 0 0 1],[1 1 1 1],[0 0 maxh maxh],[0.85 0.85 0.85])
% alpha(0.15), 
hold on
f(1) = scatter3(vpx(1,:),vpy(1,:),vpz(1,:),400,'^');
set(f(1),'MarkerEdgeColor','k',...
        'MarkerFaceColor',[1 0.9 0])
f(2) = scatter3(vpx(2,:),vpy(2,:),vpz(2,:),400,'^');
set(f(2),'MarkerEdgeColor','k',...
        'MarkerFaceColor',[1 0 0.8])
f(3) = scatter3(vpx(3,:),vpy(3,:),vpz(3,:),400,'^');
set(f(3),'MarkerEdgeColor','k',...
        'MarkerFaceColor',[0 0.8 1])
for i = 2:n %loop through layers in the structure DATA
    nn = length(DATA.(fn{i}).Rcenter);
    for j = 1:nn %loop through rows for the Radius and Center points of each layer
        x = DATA.(fn{i}).Rcenter(j,1);
        y = DATA.(fn{i}).Rcenter(j,2);
        r = mean(DATA.(fn{i}).Rradii(j,:));
        %creates a circle: increment of 0.01 is for display purposes only
        ang=0:0.01:2*pi;
        xp=r*cos(ang);
        yp=r*sin(ang);
        Hm = repmat(H(i),1,length(xp));
        p(j) = plot3(x+xp,y+yp,Hm);
        %determine the color of the circle by its radius
        [~,id] = min(abs(diams-r));
        set(p(j),'Color',c(id,:))
    end
end
axis([0 1 0 1 0 maxh])
caxis([str2double(sprintf('%.3f',mind)) str2double(sprintf('%.3f',maxd))])
leg = legend([f(3) f(2) f(1)],'Day 3','Day 2','Day 1');
set(gca,'position',[0.1 0.1 0.7 0.8])
set(leg,'box','off','position',[0.68 0.75 0.025 0.025],'FontSize',18)
% zlabel('Height Above Bed (m)','FontSize',16)
% ylabel('Along-shore Distance (m)','FontSize',16)
cb = colorbar('eastoutside');
set(cb,'position',[0.825 0.1 0.02 0.8],'LineWidth',1.5)
ylabel(cb,'Pneumatophore Diameter (m)','FontSize',16)
set(gca,'LineWidth',1.5,'FontSize',16)
grid on
hold off
view(0,90)
% export_fig([fdir 'Q2B_instLoc_oblique'],'-pdf','-nocrop','-m1','-r900')