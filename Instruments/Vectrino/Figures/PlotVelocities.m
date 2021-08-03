%Plot xyz velocities of VecPro data for an example figure for the ONR Talk.
clear
savefigs = 1;
dname = 'FSS32'; %deployment
name = 'ADHR3'; %instrument
savefigdir = 'c:\Users\bkn5\Projects\Documents\ONRTalk\';
load VP2_080315.mat
[VPRO.Data,VPRO.Config] = VPro_coordinateTransform(VPRO.Data,VPRO.Config,'bx');
start = datenum(2015,03,08,10,15,00);
stop = datenum(2015,03,08,10,15,30);
ind = find(VPRO.Data.Time >= start & VPRO.Data.Time <= stop);
time = VPRO.Data.Time(ind);
x = VPRO.Data.Profiles_VelX(ind,:);
y = VPRO.Data.Profiles_VelY(ind,:);
z1 = VPRO.Data.Profiles_VelZ1(ind,:);
z2 = VPRO.Data.Profiles_VelZ2(ind,:);
z = (z1+z2)/2;
rb = VPRO.Data.Profiles_Range;
sec = datenum(0,0,0,0,0,1);
clear VPRO

f1 = figure;
colormap jet
set(f1,'PaperOrientation','portrait',...
    'position',[400 200   1000   800]);
set(gcf,'color','w','PaperPositionMode','auto')
ax(1) = subplot(3,1,1);
imagesc(time,rb,x')
cb(1) = colorbar;
caxis([-0.4 0.4])

ax(2) = subplot(3,1,2);
imagesc(time,rb,y')
cb(2) = colorbar;
ylabel('\bf\itRange (m)','FontSize',14)
caxis([-0.4 0.4])

ax(3) = subplot(3,1,3);
imagesc(time,rb,z')
cb(3) = colorbar;
caxis([-0.1 0.1])

datetick('x','ss','keepticks','keeplimits')
xlabel('\bf\itTime (s)','FontSize',14)
set(ax(1),'position',[0.1 0.68 0.79 0.25])
set(ax(2),'position',[0.1 0.39 0.79 0.25])
set(ax(3),'position',[0.1 0.1 0.79 0.25])

set([ax(1) ax(2)],'XTickLabel',[])
set(ax,'YTick',[0.04:0.01:0.07],...
    'XLim',[start stop],'XTick',start:(sec*5):stop,...
    'FontSize',14,'FontName','Helvetica')
if savefigs
    fpath = savefigdir;fname = [dname '_' name '_rawvels'];
    export_fig([fpath fname],'-png','-native')
    disp(['Figure ' fname '.pdf saved'])
end