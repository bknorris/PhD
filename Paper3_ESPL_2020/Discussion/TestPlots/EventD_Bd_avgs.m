clear, close all
load('e:\GradSchool\DataAnalysis\Paper3\Wavelet\05-03-15\waveIGevents.mat')
sfdir = 'e:\GradSchool\DataAnalysis\Paper3\WorkingFigures\Wavelet\';
expname = 'F2F2';
exploc = {'Mudflat';'Fringe';'Forest'};
p1 = zeros(3,1);
p2 = zeros(3,1);
ltxt = cell(3,1);
fn = fieldnames(eventd);
cl = [207 176 126;
    60 166 74;
    4 76 41]./255;
cc = [180 180 180;124 124 124;51 51 51]./255;
symb = {'o';'^';'s'};
f1 = figure(1);
    set(f1,'PaperOrientation','portrait',...
        'position',[400 100   800   400],...
        'renderer','painters')
for i = 1:3
    %normalize event length by longest event, bin median bed level by event
    %length
    eventl = [eventd.(fn{i}).wave.eventl; eventd.(fn{i}).ig.eventl];
    eventl = eventl(~isnan(eventl));
    nevtl = linspace(min(eventl./max(eventl)),max(eventl./max(eventl)),30);
    bdmed = [eventd.(fn{i}).wave.deltbd; eventd.(fn{i}).ig.deltbd]; 
    bdmed = bdmed(~isnan(bdmed));
%     bdl = bdmed./eventd.(fn{i}).wave.bdst;
    [b,n,s] = bindata(eventl./max(eventl),bdmed,nevtl);
    p1(i) = errorbar(nevtl,b,s,symb{i},...
        'color',cl(i,:),...
        'linewidth',1.5,...
        'markersize',8,...
        'markeredgecolor','k',...
        'markerfacecolor',cl(i,:));hold on
    %best-fit line
    [ix,iy] = find(~isnan(b));
    pf = polyfit(nevtl(ix)',b(ix),1);
    pv = polyval(pf,nevtl(ix));
    p2(i) = plot(nevtl(ix),pv,'-',...
        'color',cc(i,:),...
        'linewidth',1.5);
    %r-squared
    yresid = b(ix)-pv';
    SSresid = sum(yresid.^2);
    SStotal = (length(b(ix))-1) * var(b(ix));
    rsq = 1 - SSresid/SStotal;
    if pf(2) < 0
        tx = '';
    else
        tx = '+';
    end
    ltxt{i} = ['y = ' sprintf('%0.1e',pf(1)) 'x' tx sprintf('%0.1e',pf(2)) ', r-sq: ' sprintf('%0.1f',rsq)];
end
txt = [exploc;ltxt];
leg = legend([p1;p2],txt{:});
%Global plot adjustments
set(gca,'ylim',[-1E-2 5E-3],'xlim',[-0.05 1.05])
ylabel('Net Bed Elevation Change (m)')
xlabel('Normalized Event Length')
prettyfigures('text',12,'labels',13,'box',1)
set(leg,'position',[0.35,0.24,0.1 0.2],...
    'fontsize',11,...
    'box','off')
% export_fig([sfdir expname '_evlBDavg'],'-png')

