%Animate Spectral fits from SlopeFromSpec
datdir = 'c:\Users\bkn5\Projects\Mekong_W2015\DataAnalysis\Spectra\Paper1\QCd\';
figdir = 'c:\Users\bkn5\Projects\Mekong_W2015\Figures\Spectra&Waves\VerticalSpectra\';
load([datdir 'HTASlopeSpectra.mat'])

%plot spectral fits, make an animation
f3 = figure(3);
set(f3,'PaperOrientation','landscape',...
    'position',[400 100   1200 600]);
set(gcf,'color','w','PaperPositionMode','auto')
%plot settings, select 'day1', 'day2', or 'day3'
dd = 'day2';%'day3';%'day1'
bins = [4 6]; %bins to plot: 1 = 1, 2 = 5, 3 = 10, 4 = 15, etc.
binname = {'bin1';'bin5';'bin10';'bin15';'bin20';'bin25';'bin35'};
name = ['HTA' dd 'Spec.gif'];
c = paruly(length(binname));
fn = fieldnames(Spec.(dd));
for j = 1:length(Spec.(dd).(fn{1}).yint)
    s(1) = subplot(131);
    for k = 1:length(bins)
        yint = Spec.(dd).(fn{1}).yint(j,bins(k));
        b = Spec.(dd).(fn{1}).sl(j,bins(k));
        f = Spec.(dd).(fn{1}).sp.f(j,:);
        sp = Spec.(dd).(fn{1}).sp.(binname{bins(k)})(j,:);
        y_hat=exp(b*log(f)+yint);
        q(k) = loglog(f,y_hat,'--k','LineWidth',1.5); hold on
        p(k) = loglog(f,sp,'+','Color',c(bins(k),:));
    end
    set(gca,'XLim',[1E0 1E1],'YLim',[1E-7 1E-3],'FontSize',12,'LineWidth',1.5)
    %         xlabel('Frequency (Hz)','FontSize',12)
    ylabel('S_z_z (m^-^2s^-^1)','FontSize',12)
    title('VP1','FontSize',12)
    hold off
    
    s(2) = subplot(132);
    for k = 1:length(bins)
        yint = Spec.(dd).(fn{2}).yint(j,bins(k));
        b = Spec.(dd).(fn{2}).sl(j,bins(k));
        f = Spec.(dd).(fn{2}).sp.f(j,:);
        sp = Spec.(dd).(fn{2}).sp.(binname{bins(k)})(j,:);
        y_hat=exp(b*log(f)+yint);
        q(k) = loglog(f,y_hat,'--k','LineWidth',1.5); hold on
        p(k) = loglog(f,sp,'+','Color',c(bins(k),:));
    end
    set(gca,'XLim',[1E0 1E1],'YLim',[1E-7 1E-3],'YTickLabel',[],'FontSize',12,'LineWidth',1.5)
    xlabel('Frequency (Hz)','FontSize',12)
    %         ylabel('S_z_z (m^-^2s^-^1)','FontSize',12)
    title('VP2','FontSize',12)
    hold off
    
    clear p
    s(3) = subplot(133);
    for k = 1:length(bins)
        yint = Spec.(dd).(fn{3}).yint(j,bins(k));
        b = Spec.(dd).(fn{3}).sl(j,bins(k));
        f = Spec.(dd).(fn{3}).sp.f(j,:);
        sp = Spec.(dd).(fn{3}).sp.(binname{bins(k)})(j,:);
        y_hat=exp(b*log(f)+yint);
        q(k) = loglog(f,y_hat,'--k','LineWidth',1.5); hold on
        p(k) = loglog(f,sp,'+','Color',c(bins(k),:),'DisplayName',binname{bins(k)});
    end
    set(gca,'XLim',[1E0 1E1],'YLim',[1E-7 1E-3],'YTickLabel',[],'FontSize',12,'LineWidth',1.5)
    
    %         xlabel('Frequency (Hz)','FontSize',12)
    %         ylabel('S_z_z (m^-^2s^-^1)','FontSize',12)
    title('VP3','FontSize',12)
    hold off
    
    suptitle(['2-Minute Spectra, ' datestr(Spec.(dd).(fn{1}).time(j),'dd/mm HH:MM')])
    set(s(1),'position',[0.08 0.12 0.28 0.75])
    set(s(2),'position',[0.38 0.12 0.28 0.75])
    set(s(3),'position',[0.68 0.12 0.28 0.75])
    leg = legend(p);
    set(leg,'LineWidth',1.5,'Position',[0.67 0.22 0.1 0.1],'box','off')
    
    drawnow
    frame = getframe(f3);
    im = frame2im(frame);
    [imind,cm] = rgb2ind(im,256);
    if j == 1;
        imwrite(imind,cm,[figdir name],'gif', 'Loopcount',inf,'DelayTime',0.5);
    else
        imwrite(imind,cm,[figdir name],'gif','WriteMode','append','DelayTime',0.5);
    end
end