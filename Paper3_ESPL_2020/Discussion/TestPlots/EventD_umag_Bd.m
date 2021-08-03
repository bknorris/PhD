%Plot event duration vs. bed movement statistics
clear, close all
load('e:\GradSchool\DataAnalysis\Paper3\Wavelet\05-03-15\waveIGevents.mat')
sfdir = 'e:\GradSchool\DataAnalysis\Paper3\WorkingFigures\Wavelet\';
expname = 'F2F2';
fn = fieldnames(eventd);
p = zeros(3,1);
for i = 1:3
    %plot event length, umag and color by delta bd
    range = 1.5E-4;
    f1 = figure(i);
    set(f1,'PaperOrientation','portrait',...
        'position',[400 100   800   400],...
        'renderer','painters');
    sp(1) = subplot(121);
    eventl = eventd.(fn{i}).wave.eventl;
    deltabd = eventd.(fn{i}).wave.deltbd;
    umag = eventd.(fn{i}).wave.umag;
    %Now find non-nan values, save row col indexes
    umag = umag(~isnan(umag));
    deltabd = deltabd(~isnan(deltabd));
    iy = linspace(-range,range,50);
    eventl = eventl(~isnan(eventl));
    cl = brewermap(length(iy),'RdBu');
    colormap(cl)
    for j = 1:length(eventl)
        [~,idx] = min(abs(iy-deltabd(j)));
        plot(eventl(j),umag(j),'o',...
            'markeredgecolor','k',...
            'markerfacecolor',cl(idx,:),...
            'markersize',6), hold on
    end
    cb(1) = colorbar;
    caxis([-range range])
    %%%
    sp(2) = subplot(122);
    eventl = eventd.(fn{i}).ig.eventl;
    deltabd = eventd.(fn{i}).ig.deltbd;
    umag = eventd.(fn{i}).ig.umag;
    %Now find non-nan values, save row col indexes
    umag = umag(~isnan(umag));
    deltabd = deltabd(~isnan(deltabd));
    iy = linspace(-2E-4,2E-4,50);
    eventl = eventl(~isnan(eventl));
    cl = brewermap(length(iy),'RdBu');
    colormap(cl)
    for j = 1:length(eventl)
        [~,idx] = min(abs(iy-deltabd(j)));
        plot(eventl(j),umag(j),'o',...
            'markeredgecolor','k',...
            'markerfacecolor',cl(idx,:),...
            'markersize',6), hold on
    end
    cb(2) = colorbar;
    caxis([-range range])
    %%%
    %Global plot adjustments
    title(sp(1),'Wave')
    title(sp(2),'IG')
    xlabel(sp(1),'Event Length (min)')
    xlabel(sp(2),'Event Length (min)')
    ylabel(sp(1),'Velocity Magnitude (m/s)')
    ylabel(sp(2),'Velocity Magnitude (m/s)')
    ylabel(cb(1),['\Delta ' 'BLE (m)'])
    ylabel(cb(2),['\Delta ' 'BLE (m)'])
    set(sp(1),'xlim',[10^-2 1],'ylim',[0 0.2],...
        'position',[0.10 0.16 0.23 0.75],...
        'xscale','log')
    set(sp(2),'xlim',[10^-2 100],'ylim',[0 0.1],...
        'position',[0.61 0.16 0.23 0.75],...
        'xscale','log')
    prettyfigures('text',12,'labels',13,'box',1)
    export_fig([sfdir expname '_' fn{i} '_evlUmagBd'],'-png')
end

