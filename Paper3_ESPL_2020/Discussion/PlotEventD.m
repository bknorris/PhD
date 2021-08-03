%Plot event duration vs. bed movement statistics
clear, close all
load('d:\Mekong_W2015\DataAnalysis\Paper3\Wavelet\05-03-15\waveIGevents.mat')
fn = fieldnames(eventd);
f1 = figure(1);
set(f1,'PaperOrientation','portrait',...
    'position',[400 100   500   400],...
    'renderer','painters');
cl = [207 176 126;
    60 166 74;
    4 76 41]./255;
cc = [180 180 180;124 124 124;51 51 51]./255;
p = zeros(3,1);
for r = 1:3
    ord = [3 2 1];
    i = ord(r);
    %plot wave
    eventl = eventd.(fn{i}).wave.eventl;
    dbd = eventd.(fn{i}).wave.deltbd;
    sp(1) = subplot(121);hold on
    wdb1 = dbd(dbd>0);
    wel1 = eventl(dbd>0);
    dd1 = max(wdb1); %scaling factor
    plot(wel1,wdb1./dd1,'o',...
        'markeredgecolor',cc(i,:),...
        'markersize',4)
    
    %plot IG
    eventl = eventd.(fn{i}).ig.eventl;
    dbd = eventd.(fn{i}).ig.deltbd;
    wdb2 = dbd(dbd>0);
    wel2 = eventl(dbd>0);
    dd2 = max(wdb2); %scaling factor
    plot(wel2,wdb2./dd2,'^',...
        'markeredgecolor',cc(i,:),...
        'markersize',6)
    
    sp(2) = subplot(122);hold on
    %plot wave
    eventl = eventd.(fn{i}).wave.eventl;
    dbd = eventd.(fn{i}).wave.deltbd;
    wdb1 = dbd(dbd<0);
    wel1 = eventl(dbd<0);
    dd1 = min(wdb1); %scaling factor
    plot(wel1,wdb1./dd1,'o',...
        'markeredgecolor',cc(i,:),...
        'markersize',4)
    pct = nnz(dbd>0)/nnz(dbd<0);
    disp(['# of accretionary to erosive events, wave: ' num2str(pct) ', ' fn{i}])

    %plot IG
    eventl = eventd.(fn{i}).ig.eventl;
    dbd = eventd.(fn{i}).ig.deltbd;
    wdb2 = dbd(dbd<0);
    wel2 = eventl(dbd<0);
    dd2 = min(wdb2); %scaling factor
    plot(wel2,wdb2./dd2,'^',...
        'markeredgecolor',cc(i,:),...
        'markersize',6)
    pct = nnz(dbd>0)/nnz(dbd<0);
    disp(['# of accretionary to erosive events, IG: ' num2str(pct) ', ' fn{i}])
end
for i = 1:3
    sp(1) = subplot(121);hold on
    %Wave
    eventl = eventd.(fn{i}).wave.eventl;
    dbd = eventd.(fn{i}).wave.deltbd;
    wdb1 = dbd(dbd>0);
    wel1 = eventl(dbd>0);
    dd1 = max(wdb1); %scaling factor
    %IG
    eventl = eventd.(fn{i}).ig.eventl;
    dbd = eventd.(fn{i}).ig.deltbd;
    wdb2 = dbd(dbd>0);
    wel2 = eventl(dbd>0);
    dd2 = max(wdb2); %scaling factor
    evl = [wel1; wel2];
    ddb = [wdb1./dd1; wdb2./dd2];
    pf = polyfit(log10(evl),log10(ddb),1);
    xs = linspace(min(evl),max(evl),100);
    ys = xs.^pf(1).*exp(pf(2));
    plot(xs,ys,'-',...
        'linewidth',2,...
        'color',cl(i,:))
    
    sp(2) = subplot(122);hold on
    %Wave
    eventl = eventd.(fn{i}).wave.eventl;
    dbd = eventd.(fn{i}).wave.deltbd;
    wdb1 = dbd(dbd<0);
    wel1 = eventl(dbd<0);
    dd1 = min(wdb1); %scaling factor
    %IG
    eventl = eventd.(fn{i}).ig.eventl;
    dbd = eventd.(fn{i}).ig.deltbd;
    wdb2 = dbd(dbd<0);
    wel2 = eventl(dbd<0);
    dd2 = min(wdb2); %scaling factor
    evl = [wel1; wel2];
    ddb = [wdb1./dd1; wdb2./dd2];
    pf = polyfit(log10(evl),log10(ddb),1);
    xs = linspace(min(evl),max(evl),100);
    ys = xs.^pf(1).*exp(pf(2));
    p(i) = plot(xs,ys,'-',...
        'linewidth',2,...
        'color',cl(i,:));
end
title(sp(1),'Accretion')
title(sp(2),'Erosion')
xlabel(sp(1),'Event Length (min)')
xlabel(sp(2),'Event Length (min)')
ylabel(sp(1),'Normalized Net Elevation Change')
set(sp(1),'xscale','log','yscale','log',...
    'xlim',[10^-2 10^2],'ylim',[10^-6 2],...
    'position',[0.13 0.14 0.38 0.79])
set(sp(2),'xscale','log','yscale','log',...
    'xlim',[10^-2 10^2],'ylim',[10^-6 2],...
    'yticklabel',[],...
    'position',[0.56 0.14 0.38 0.79])
leg = legend(p,{'Mudflat','Fringe','Forest'});
prettyfigures('text',12,'labels',13,'box',1)
set(leg,'position',[0.79 0.22 0.05 0.05])
sfdir = 'd:\Projects\Mekong_W2015\Figures\Paper3\Discussion\';
% export_fig([sfdir 'F2F2_Acc_Ero_nbd_evtl'],'-png')

f2 = figure(2);
set(f2,'PaperOrientation','portrait',...
    'position',[400 100   500   400],...
    'renderer','painters');hold on
cl = [207 176 126;
    60 166 74;
    4 76 41]./255;
p = zeros(3,1);
for i = 1:3
    %plot wave
    eventl = eventd.(fn{i}).wave.eventl;
    bdvar = eventd.(fn{i}).wave.bdvar;
    plot(eventl,bdvar,'*',...
        'markeredgecolor',cl(i,:),...
        'markersize',6)
    %plot IG
    eventl = eventd.(fn{i}).ig.eventl;
    bdvar = eventd.(fn{i}).ig.bdvar;
    z = plot(eventl,bdvar,'^',...
        'markerfacecolor',cl(i,:),...
        'markersize',6,...
        'markeredgecolor','k');
    p(i) = z(1);
end
set(gca,'xscale','log','yscale','log',...
    'xlim',[10^-3 10^2],...
    'ylim',[10^-14 10^-2])
leg = legend(p,{'Mudflat','Fringe','Forest'});
xlabel('Event Length (min)')
ylabel('\sigma^2_{bd}/\Deltat*U')
sfdir = 'd:\Projects\Mekong_W2015\Figures\Paper3\Discussion\';
prettyfigures('text',12,'labels',13,'box',1)
set(leg,'position',[0.75 0.22 0.05 0.05])
% export_fig([sfdir 'F2F2_Bdnorm_evtl'],'-png')
