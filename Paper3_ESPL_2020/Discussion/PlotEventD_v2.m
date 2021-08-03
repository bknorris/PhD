%Plot event duration vs. bed movement statistics
clear, close all
load('e:\GradSchool\DataAnalysis\Paper3\Wavelet\05-03-15\waveIGevents.mat')
fn = fieldnames(eventd);
%first plot the # of wave vs. IG events by period!
cl = [207 176 126;
    60 166 74;
    4 76 41]./255;
cc = [180 180 180;124 124 124;51 51 51]./255;
p = zeros(3,1);
for i = 1:3
    f1 = figure(i);
    set(f1,'PaperOrientation','portrait',...
        'position',[400 100   500   400],...
        'renderer','painters');
    eventl = eventd.(fn{i}).wave.eventl;
    deltabd = eventd.(fn{i}).wave.deltbd;
    umag = eventd.(fn{i}).wave.umag;
    period = eventd.(fn{i}).wave.period;
    
    
    
    
    
    
    
end
title(sp(1),'Accretion')
title(sp(2),'Erosion')
xlabel(sp(1),'Event Length (min)')
xlabel(sp(2),'Event Length (min)')
ylabel(sp(1),'Normalized Net Elevation Change')
set(sp(1),'xlim',[-0.09 10],'ylim',[0 10^-3],...
    'position',[0.13 0.14 0.38 0.79])
set(sp(2),'xlim',[-0.09 5],'ylim',[-10^-4 0],...
    'yticklabel',[],...
    'position',[0.56 0.14 0.38 0.79])
leg = legend(p,{'Mudflat','Fringe','Forest'});
prettyfigures('text',12,'labels',13,'box',1)
set(leg,'position',[0.79 0.22 0.05 0.05])
sfdir = 'd:\Projects\Mekong_W2015\Figures\Paper3\Discussion\';
% export_fig([sfdir 'F2F2_Acc_Ero_nbd_evtl'],'-png')
