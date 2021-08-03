%Plot event duration vs. bed movement statistics
clear, close all
load('d:\Mekong_W2015\DataAnalysis\Paper3\Wavelet\05-03-15\waveIGevents.mat')
sfdir = 'e:\GradSchool\DataAnalysis\Paper3\WorkingFigures\Wavelet\';
expname = 'F2F2';    
% f1 = figure(1);
%     set(f1,'PaperOrientation','portrait',...
%         'position',[400 100   800   400],...
%         'renderer','painters');
fn = fieldnames(eventd);
p = zeros(3,1);
for i = 1:3
    figure
%     sp(i) = subplot(1,3,i);
    %plot event length, deltabd and color by period
    eventl = [eventd.(fn{i}).wave.eventl; eventd.(fn{i}).ig.eventl];
    phase = [eventd.(fn{i}).wave.phase; eventd.(fn{i}).ig.phase];
    period = [eventd.(fn{i}).wave.period eventd.(fn{i}).ig.period];
    %Now find non-nan values, save row col indexes
    [ix,~] = find(~isnan(phase));
    phase = phase(~isnan(phase));
    idx = linspace(min(phase),max(phase),length(phase));
    eventl = eventl(~isnan(eventl));
    nevntl = eventl./max(eventl);
    cl = brewermap(length(idx),'RdBu');
    for j = 1:length(eventl)
        [~,iy] = min(abs(idx-phase(j)));
        plot(nevntl(j),period(ix(j)),'o',...
            'markeredgecolor','k',...
            'markerfacecolor',cl(iy,:),...
            'markersize',6), hold on
    end 
    colormap(cl)
    cb(1) = colorbar;
    caxis([min(idx) max(idx)])
    set(gca,'xscale','log','yscale','log')
end

