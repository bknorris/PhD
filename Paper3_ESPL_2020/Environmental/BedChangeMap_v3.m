clear, close all
%plot bed map
figdir = 'd:\GradSchool\Writing\GM_2018_Norris\Figures\V5\';
f1 = figure(1);
set(f1,'PaperOrientation','portrait',...
    'renderer','painters',...
    'position',[400 100   450   900]);
dirc = 'e:\Mekong_W2015\Images\Quadrats\Q2\Reconstructions\';
qdat = 'e:\Mekong_W2015\Images\Quadrats\Q2\Reconstructions\2A\Analysis\';
fname = {'twobraster.txt','twobtwocdiff.txt','twocraster.txt'};
%Load bed change data
vpxy = [0.8033 0.245;0.5997 0.245;0.4978 0.245];
rotx = [0,0];
roty = [0,0];
rotz = [0,0];
cc = brewermap(15,'Greys');symb = {'o';'d';'p'};
cl = flipud([0.1 0.1 0.1;0.5 0.5 0.5;0.7 0.7 0.7]);
colormap(cc);
for i = 1:3
    [z,r] = ascread([dirc fname{i}]);
    z(z<-9000) = NaN;
    %[z,r] = arcgridread([dirc fname{i}]);
    ccscale = 6.657;
    Z = (z./ccscale).*1000; %*1000 converts to mm
    if i == 1
        Z = Z+895.6845; %adjust for difference between lowest point of twoCraster and twoBraster
    elseif i == 3
        Z = Z+892.56; %adjust for distance under VP2 based on end of bdtrace timeseries
    end
    [m,n] = size(Z);
    x = linspace(0,1,n);X = repmat(x,m,1);
    y = linspace(1,0,m);Y = repmat(y,n,1)';
    %load quadrat data
    load([qdat 'Q2Astats_fixed.mat'])
    
    %contour using contourf
    sp(i) = subplot(3,1,i);
    [c,h] = contourf(X,Y,Z);
    set(h,'EdgeColor','none')
    hold on
    
    %calculate radii of pneumatophores at level 1 of the reconstruction
    nn = length(DATA.layer1.Rcenter);
    for j = 1:nn %loop through rows for the Radius and Center points of each layer
        cx = DATA.layer1.Rcenter(j,1);
        cy = DATA.layer1.Rcenter(j,2);
        %transformation matrix, rotate 90deg CCW
        R = [0 -1;1 0];
        rot = R*[cx;cy];
        cx = abs(rot(1));cy = rot(2);
        cr = mean(DATA.layer1.Rradii(j,:));
        ang=0:0.01:2*pi;
        xp=cr.*cos(ang);
        yp=cr.*sin(ang);
        %plot circles with the same radius as the vegetation
        H = patch(cx+xp,cy+yp,1);
        set(H,'FaceColor',[0 0 0],'EdgeColor','none')
    end
    for j = 1:3
        plot(vpxy(j,1),vpxy(j,2),symb{j},...
            'Color','k',...
            'MarkerSize',10,...
            'MarkerFaceColor',cl(j,:),...
            'LineWidth',1.5)
    end
    if i == 2
        caxis([-30 30]),
    else
        caxis([-50 50])
    end
    grid on
    set(gca,'gridlinestyle',':')
end
%labels
cb1 = colorbar('peer',sp(1),'location','eastoutside');
cb2 = colorbar('peer',sp(2),'location','eastoutside');
cb3 = colorbar('peer',sp(3),'location','eastoutside');
ylabel(cb1,'Bed Level [mm]')
ylabel(cb2,'Bed Level Change [mm]')
ylabel(cb3,'Bed Level [mm]')
ylabel(sp(1),'Along-shore Distance [m]')
ylabel(sp(2),'Along-shore Distance [m]')
ylabel(sp(3),'Along-shore Distance [m]')
xlabel(sp(3),'Across-shore  Distance [m]')

%positioning
set(sp(1),'position',[0.16 0.73 0.6 0.25],...
    'xtick',0:0.25:1,'ytick',0:0.25:1)
set(sp(2),'position',[0.16 0.405 0.6 0.25],...
    'xtick',0:0.25:1,'ytick',0:0.25:1)
set(sp(3),'position',[0.16 0.08 0.6 0.25],...
    'xtick',0:0.25:1,'ytick',0:0.25:1)
set(cb1,'position',[0.81 0.73 0.05 0.25],'ytick',[-50 -25 0 25 50])
set(cb2,'position',[0.81 0.405 0.05 0.25],'ytick',[-30 -15 0 15 30])
set(cb3,'position',[0.81 0.08 0.05 0.25],'ytick',[-50 -25 0 25 50])
prettyfigures('text',13,'labels',14,'box',1)
set(f1,'units','inches');
pos = get(f1,'Position');

set(f1,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(f1,[figdir 'BedChangeMap_v3grey'],'-dpdf','-r0')