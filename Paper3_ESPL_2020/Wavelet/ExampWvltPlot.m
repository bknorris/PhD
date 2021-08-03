%Plot example wavelets for the manuscript
clear,close all
dir1 = 'e:\Mekong_W2015\DataAnalysis\Paper3\Wavelet\05-03-15\';
files = dir([dir1 '*wvlt.mat']);files = {files.name};
f1 = figure(1);
set(f1,'PaperOrientation','portrait',...
    'position',[400 100   1000   600],...
    'renderer','painters');
cb = brewermap(100,'*RdYlBu');
sp = zeros(4,1);
llim = -5.86;ulim = -4;
for i = 1:3
    load([dir1 files{i}])
    dat = wvlt.x;clear wvlt
    %Plot Routine
    sp(i) = subplot(1,4,i);
    plotwtc(dat.Rsq,dat.period,dat.coi,dat.sig95,dat.t,dat.Wxy,dat.dt)
        colormap(gca,cb)
    if i == 2
        hold on,
        %plot rectangle showing inset
        x1 = 45;x2 = 50;
        y1 = llim;y2 = ulim;
        x = [x1, x2, x2, x1, x1];
        y = [y1, y1, y2, y2, y1];
        plot(x, y,'-','color','r', 'LineWidth', 2);
        %inset window
        sp(4) = subplot(144); 
        plotwtc(dat.Rsq,dat.period,dat.coi,dat.sig95,dat.t,dat.Wxy,dat.dt)
        colormap(gca,cb)
        ax = gca;
    end
    clear dat
end
%Plot Adjustments
set(sp(1),'position',[0.11 0.23 0.26 0.4])
set(sp(2),'position',[0.405 0.23 0.26 0.4],...
    'yticklabel',[])
set(sp(3),'position',[0.7 0.23 0.26 0.4],...
    'yticklabel',[])
%Inset plot
Yticks = 2.^(fix(llim):fix(ulim));
set(sp(4),'position',[0.48 0.73 0.4 0.25],...
    'ylim',[-5.86 -4],'xlim',[45 50],...
    'YTick',log2(Yticks(:)), ...
    'YTickLabel',num2str(Yticks'.*60))
set([sp(1) sp(2) sp(3)],'xlim',[0 80],'xtick',0:20:80)
set(ax,{'XColor','YColor'},{'r','r'})
%Labels
xlabel(sp(1),'Minutes elapsed after 15:30')
xlabel(sp(2),'Minutes elapsed after 15:30')
xlabel(sp(3),'Minutes elapsed after 15:30')
% suplabel('Time on 05-03-15','x')
ylabel(sp(1),'Period [min]')
xlabel(sp(4),'Minutes elapsed after 15:30')
ylabel(sp(4),'Period [sec]')
cbar = colorbar('location','southoutside');
set(cbar,'position',[0.24 0.1 0.6 0.03],...
    'linewidth',1.5,'tickdir','out')
xlabel(cbar,'Magnitude Squared Coherence','fontsize',13)
prettyfigures('text',12,'labels',13,'box',1)
sfdir = 'd:\GradSchool\DataAnalysis\Paper3\Figures\';
set(f1,'units','inches');
pos = get(f1,'Position');
% set(f1,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
% print(f1,[sfdir 'WvltExample'],'-dpdf','-r0')