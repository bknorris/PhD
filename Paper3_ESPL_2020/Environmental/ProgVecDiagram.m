%Progressive sediment flux vectors for Paper #3
clear, close all
dir1 = 'd:\Mekong_W2015\DataAnalysis\TOS\';
UWaqd = {'AQD_SW_2014.mat';'AQD_SW_2015.mat'};
scale = [2.5E-6; 1.2E-6];
lat = [9.494083 9.4930472];
lon = [106.2435 106.243608];
%Plot Routine
f1 = figure(1);
set(f1,'PaperOrientation','portrait',...
    'position',[400 100 600 400]);
set(gcf,'color','w','PaperPositionMode','auto','Renderer','OpenGL')
hold on
fdir = 'd:\Mekong_W2015\DataAnalysis\Paper3\';
load([fdir 'CLD_landboundaries.mat']);
%Plot the landboundaries file
for i = 6:11
    x = dat.forest(i).X;
    y = dat.forest(i).Y;
    plot(x,y,'k','linewidth',1.5)
end
x = dat.fringe.X;
y = dat.fringe.Y;
plot(x,y,'k','linewidth',1.5)
for i = 2
    x = dat.urban(i).X;
    y = dat.urban(i).Y;
    plot(x,y,'k','linewidth',1.5)
end
%Plot progressive vectors
for ff = 1:length(UWaqd)
    load([dir1 UWaqd{ff}])
    %remove out of air values
    aqd.v1(aqd.air) = NaN;
    aqd.v2(aqd.air) = NaN;
    u = nanmean(aqd.v1,2);
    v = nanmean(aqd.v2,2);
    %calc SSC (from Aaron)
    SSC = mean([0.0688.*aqd.ext1+50 0.0688.*aqd.ext2+50],2);
    uc = u.*SSC;vc = v.*SSC;
    innX = isnan(uc);
    innY = isnan(vc);
    uc(innX) = 0;
    vc(innY) = 0;
    c = brewermap(length(uc),'Reds');
    for i = 1:length(uc)
        umean(i) = uc(i);
        vmean(i) = vc(i);
        posX = cumsum([lon(ff) ; umean(:).*scale(ff)]);
        posY = cumsum([lat(ff) ; vmean(:).*scale(ff)]);
    end
    posX([isnan(1) ; innX(:)]) = NaN;  % NaNs re-inserted in their original locations.
    posY([isnan(1) ; innY(:)]) = NaN;
    for i = 1:length(uc)-1
        line(posX(i:i+1),posY(i:i+1),...
            'LineWidth',1.5,...
            'Color',c(i,:));
    end
end
%colorbar properties
colormap(c);
cb = colorbar;
caxis([1 14])
ylabel(cb,'Days of Experiments','fontsize',13)
set(cb,'linewidth',1.5,'tickdir','out')
%global adjustments
set(gca,'xlim',[106.16 106.29],...
    'ylim',[9.46 9.56])
xlabel('Longitude (dd)')
ylabel('Latitude (dd)')
prettyfigures('text',12,'labels',13,'box',1)
figdir = 'e:\GradSchool\DataAnalysis\Paper3\Figures\';
export_fig([figdir 'ProgVecDiag_v2'],'-pdf','-nocrop')

