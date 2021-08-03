clear, close all
load('D:\Projects\Mekong_W2015\DataAnalysis\Paper2\RS_estimates\ExampleRS_HTA1vp3_v2.mat')
fdir = 'g:\GradSchool\DataAnalysis\Paper2\WorkingFigures\ReynoldsStress\CF_V4\';

%A good fit: Use i = 1, ii = 3, j = 3, jj = 8;
f1 = figure(1);
set(f1,'PaperOrientation','portrait',...
    'position',[400 100   900   400]);
set(gcf,'color','w','paperpositionmode','auto') 
dk = k(2)-k(3);
%%
sp(1) = subplot(131);
xx = [flipud(k);k];
yy = [flipud(COuwerrh.*k);COuwerrl.*k];
fill(xx(2:end-1),-1*yy(2:end-1),[0.7 0.7 0.7],...
   'EdgeColor','none'), hold on
p1(1) = plot(k,-1*COuw.*k,'color',[0.7 0.7 0.7],'linewidth',1);
p1(2) = plot(k,-COuwstar.*uw.*k,'k','linewidth',1.5);
plot(ones(10,1)*kc,linspace(-0.1,0.1,10),'color','k')
rectangle('position',[1*10^2,-5E-3,1.9E3,12E-3])
txt = 'kCo_{u''w''} (rad m/s^2)';
text(kc+20,0.02,'k_c')
xlabel('k (rad/m)')
ylabel(txt)
%%
sp(2) = subplot(132);
%bin data
bins = linspace(min(k),max(k),200);
COuwk = bindata(k,-1*COuw.*k,bins);
p1(1) = plot(bins,COuwk,'color',[0.7 0.7 0.7],'linewidth',1);hold on
p1(2) = plot(k,-COuwstar.*uw.*k,'k','linewidth',1.5);
k0 = 320;
plot(ones(10,1)*k0,linspace(-1E-3,10E-3,10),'color','k')
plot(ones(10,1)*kc,linspace(-1E-3,10E-3,10),'color','k')
text(kc+10,1.1E-3,'k_c')
text(k0+20,1.1E-3,'k_0')
%%
sp(3) = subplot(133);
p2(1) = semilogx(k,-1*cumsum(COuw.*dk),'color',[0.7 0.7 0.7],'linewidth',1.5);hold on
p2(2) = semilogx(k,-1*cumsum(COuwstar.*uw.*dk),'k','linewidth',1.5);hold on
plot(ones(10,1)*kc,linspace(-1E-3,10E-3,10),'color','k')
text(kc+90,3E-3,'k_c')
txt = '$\int{Co_{u''w''}\enspace dk}\quad (m/s)$';
xlabel('k (rad/m)')
ylabel(txt,'interpreter','latex')
%%
set(sp(1),...
    'xscale','log',...
    'xlim',[10^-0.3 10^3.3],...
    'xtick',[10^0 10^1 10^2 10^3],...
    'ylim',[-0.03 0.05],...
    'ytick',-0.02:0.02:0.04,...
    'position',[0.1 0.15 0.36 0.76])
set(sp(2),...
    'xscale','log',...
    'xlim',[10^2 10^3.3],...
    'xtick',[10^2 10^3],...
    'ylim',[-1E-3 1.5E-3],...
    'ytick',-1E-3:1E-3:1E-3,...
    'position',[0.4 0.6 0.18 0.33])
set(sp(3),...
    'xdir','reverse',...
    'xlim',[10^-0.3 10^3.3],...
    'xtick',[10^0 10^1 10^2 10^3],...
    'ylim',[-1E-4 4E-3],...
    'ytick',0:2E-3:6E-3,...
    'position',[0.64 0.15 0.35 0.76])
prettyfigures('text',12,'labels',13,'box',1)
export_fig([fdir 'VP3_HTA1_goodfit_v3'],'-pdf')