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
    eventl = eventd.(fn{i}).wave.eventl;
    bdvar = eventd.(fn{i}).wave.bdvar;
    period = eventd.(fn{i}).wave.period;
    %Now find non-nan values, save row col indexes
    [iy,~] = find(~isnan(bdvar));
    bdvar = bdvar(~isnan(bdvar));
    eventl = eventl(~isnan(eventl));
    cl = brewermap(length(period),'RdBu');
    colormap(cl)
    for j = 1:length(bdvar)
        plot(eventl(j),bdvar(j),'o',...
            'markeredgecolor','k',...
            'markerfacecolor',cl(iy(j),:),...
            'markersize',6), hold on
    end
    cb(1) = colorbar;
    caxis([min(period*60) max(period*60)])
    %%%
    sp(2) = subplot(122);
    eventl = eventd.(fn{i}).ig.eventl;
    bdvar = eventd.(fn{i}).ig.bdvar;
    umag = eventd.(fn{i}).ig.umag;
    period = eventd.(fn{i}).ig.period;
    %Now find non-nan values, save row col indexes
    [iy,~] = find(~isnan(bdvar));
    bdvar = bdvar(~isnan(bdvar));
    eventl = eventl(~isnan(eventl));
    cl = brewermap(length(period),'RdBu');
    for j = 1:length(bdvar)
        plot(eventl(j),bdvar(j),'o',...
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
    xlabel(sp(1),'Event Length (min)')
    xlabel(sp(2),'Event Length (min)')
    ylabel(sp(1),'\sigma_{b}/\Deltat*U')
    ylabel(sp(2),'\sigma_{b}/\Deltat*U')
    ylabel(cb(1),'Period (sec)')
    ylabel(cb(2),'Period (min)')
    set(sp(1),'xlim',[10^-2 1],'ylim',[10^-12 10^-2],...
        'position',[0.11 0.16 0.25 0.75],...
        'xscale','log','yscale','log')
    set(sp(2),'xlim',[10^-2 100],'ylim',[10^-12 10^-2],...
        'position',[0.6 0.16 0.25 0.75],...
        'xscale','log','yscale','log')
    prettyfigures('text',12,'labels',13,'box',1)
	export_fig([sfdir expname '_' fn{i} '_evlBDvarPrd'],'-png')
end

