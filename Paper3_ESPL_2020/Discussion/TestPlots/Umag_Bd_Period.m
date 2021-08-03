%Plot event duration vs. bed movement statistics
clear, close all
load('e:\GradSchool\DataAnalysis\Paper3\Wavelet\05-03-15\waveIGevents.mat')
sfdir = 'e:\GradSchool\DataAnalysis\Paper3\WorkingFigures\Wavelet\';
expname = 'F2F2';
fn = fieldnames(eventd);
p = zeros(3,1);
for i = 1:3
    %plot event length, deltabd and color by period
    f1 = figure(i);
    set(f1,'PaperOrientation','portrait',...
        'position',[400 100   800   400],...
        'renderer','painters');
    sp(1) = subplot(121);
    umag = eventd.(fn{i}).wave.umag;
    deltabd = eventd.(fn{i}).wave.deltbd;
    period = eventd.(fn{i}).wave.period;
    %Now find non-nan values, save row col indexes
    [iy,~] = find(~isnan(deltabd));
    deltabd = deltabd(~isnan(deltabd));
    umag = umag(~isnan(umag));
    cl = brewermap(length(period),'RdBu');
    colormap(cl)
    for j = 1:length(umag)
        plot(umag(j),deltabd(j),'o',...
            'markeredgecolor','k',...
            'markerfacecolor',cl(iy(j),:),...
            'markersize',6), hold on
    end
    cb(1) = colorbar;
    caxis([min(period*60) max(period*60)])
    %%%
    sp(2) = subplot(122);
    umag = eventd.(fn{i}).ig.umag;
    deltabd = eventd.(fn{i}).ig.deltbd;
    period = eventd.(fn{i}).ig.period;
    %Now find non-nan values, save row col indexes
    [iy,~] = find(~isnan(deltabd));
    deltabd = deltabd(~isnan(deltabd));
    umag = umag(~isnan(umag));
    cl = brewermap(length(period),'RdBu');
    colormap(cl)
    for j = 1:length(umag)
        plot(umag(j),deltabd(j),'o',...
            'markeredgecolor','k',...
            'markerfacecolor',cl(iy(j),:),...
            'markersize',6), hold on
    end
    cb(2) = colorbar;
    caxis([min(period) max(period)])
    %%%
    %Global plot adjustments
    title(sp(1),'Wave')
    title(sp(2),'IG')
    xlabel(sp(1),'Velocity Magnitude (m/s)')
    xlabel(sp(2),'Velocity Magnitude (m/s)')
    ylabel(sp(1),'BLE (m)')
    ylabel(sp(2),'BLE (m)')
    ylabel(cb(1),'Period (sec)')
    ylabel(cb(2),'Period (min)')
    set(sp(1),'xlim',[10^-2 1],'ylim',[-5E-4 5E-4],...
        'position',[0.11 0.16 0.25 0.75],...
        'xscale','log')
    set(sp(2),'xlim',[10^-2 1],'ylim',[-5E-4 5E-4],...
        'position',[0.6 0.16 0.25 0.75],...
        'xscale','log')
    prettyfigures('text',12,'labels',13,'box',1)
    export_fig([sfdir expname '_' fn{i} '_umagBdPrd'],'-png')
end

