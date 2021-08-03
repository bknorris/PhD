%Plot the combined event data
clear,close all
load('G:\Mekong_W2015\DataAnalysis\Paper3\Wavelet\CombEventData.mat');
units = {'min';'m';'m';'m/s';'m/s'};
cl = [207 176 126;
    60 166 74;
    4 76 41]./255;
symb = {'o';'^';'s'};
fn = fieldnames(dat);
%% Median Bed Level
eb = zeros(3,1);
for i = 1:length(fn)
    %Try plotting uorb, umag, eventl against the bed movement stats
    %event length
    f1 = figure(1);
    set(f1,'PaperOrientation','portrait',...
        'position',[400 100   800   400],...
        'renderer','painters')
    eventl = dat.(fn{i}).eventl;
    %     nevtl = eventl./max(eventl);
    bdmed = dat.(fn{i}).bdmed;
    bins = linspace(min(eventl),max(eventl),10);
    sp(1) = subplot(211);
    [b,~,s] = bindata(eventl,bdmed,bins);
    eb(i) = errorbar(bins,b,s,symb{i},...
        'color',cl(i,:),...
        'linewidth',1.5,...
        'markersize',8,...
        'markeredgecolor','k',...
        'markerfacecolor',cl(i,:));hold on
    sp(2) = subplot(212);
    [n1,x1] = hist(eventl,bins);
    bar(x1,n1./numel(eventl),...
        'facecolor',cl(i,:),...
        'edgecolor','none'),hold on
end
ytext = sprintf(['Median elev. relative' '\n' 'to event beginning [m]']);
xlabel(sp(2),'Event Length [min]')
ylabel(sp(2),'% occurence')
ylabel(sp(1),ytext)
leg = legend(eb,fn);
set(leg,'position',[0.8 0.32 0.05 0.05])
%%%
eb = zeros(3,1);
for i = 1:length(fn)
    %Try plotting uorb, umag, eventl against the bed movement stats
    %uorb - orbital wave velocity
    f2 = figure(2);
    set(f2,'PaperOrientation','portrait',...
        'position',[400 100   800   400],...
        'renderer','painters')
    uorb = dat.(fn{i}).uorb;
    bdmed = dat.(fn{i}).bdmed;
    bins = linspace(min(uorb),max(uorb),10);
    sp(1) = subplot(211);
    [b,~,s] = bindata(uorb,bdmed,bins);
    eb(i) = errorbar(bins,b,s,symb{i},...
        'color',cl(i,:),...
        'linewidth',1.5,...
        'markersize',8,...
        'markeredgecolor','k',...
        'markerfacecolor',cl(i,:));hold on
    sp(2) = subplot(212);
    [n1,x1] = hist(uorb,bins);
    bar(x1,n1./numel(uorb),...
        'facecolor',cl(i,:),...
        'edgecolor','none'),hold on
end
xlabel('Orbital Wave Velocity [m/s]')
ytext = sprintf(['Median elev. relative' '\n' 'to event beginning [m]']);
ylabel(sp(2),'% occurence')
ylabel(sp(1),ytext)
leg = legend(eb,fn);
set(leg,'position',[0.8 0.32 0.05 0.05])
%%%
eb = zeros(3,1);
for i = 1:length(fn)
    %Try plotting uorb, umag, eventl against the bed movement stats
    %umag - cross shore velocity magnitude
    f3 = figure(3);
    set(f3,'PaperOrientation','portrait',...
        'position',[400 100   800   400],...
        'renderer','painters')
    umag = dat.(fn{i}).umag;
    bdmed = dat.(fn{i}).bdmed;
    bins = linspace(min(umag),max(umag),10);
    sp(1) = subplot(211);
    [b,~,s] = bindata(umag,bdmed,bins);
    eb(i) = errorbar(bins,b,s,symb{i},...
        'color',cl(i,:),...
        'linewidth',1.5,...
        'markersize',8,...
        'markeredgecolor','k',...
        'markerfacecolor',cl(i,:));hold on
    sp(2) = subplot(212);
    [n1,x1] = hist(umag,bins);
    bar(x1,n1./numel(umag),...
        'facecolor',cl(i,:),...
        'edgecolor','none'),hold on
end
xlabel('Cross-shore velocity magnitude [m/s]')
ytext = sprintf(['Median elev. relative' '\n' 'to event beginning [m]']);
ylabel(sp(2),'% occurence')
ylabel(sp(1),ytext)
leg = legend(eb,fn);
set(leg,'position',[0.8 0.32 0.05 0.05])
%%%
%% Net Bed Change
% eb = zeros(3,1);
% for i = 1:length(fn)
%     %Try plotting uorb, umag, eventl against the bed movement stats
%     %event length
%     f4 = figure(4);
%     set(f4,'PaperOrientation','portrait',...
%         'position',[400 100   800   400],...
%         'renderer','painters')
%     eventl = dat.(fn{i}).eventl;
%     deltbd = dat.(fn{i}).deltbd;
%     sp(1) = subplot(211);
%     id = find(deltbd>0);
%     bins = linspace(min(eventl(id)),max(eventl(id)),10);
%     [b,~,s] = bindata(eventl,deltbd(id),bins);
%     U = s;
%     elim = find(b-s<=0);
%     L = s;
%     L(elim) = b(elim);
%     eb(i) = errorbar(bins,b,L,U,symb{i},...
%         'color',cl(i,:),...
%         'linewidth',1.5,...
%         'markersize',8,...
%         'markeredgecolor','k',...
%         'markerfacecolor',cl(i,:));hold on
%     sp(2) = subplot(212);
%     id = find(deltbd<0);
%     bins = linspace(min(eventl(id)),max(eventl(id)),10);
%     [b,~,s] = bindata(eventl,deltbd(id),bins);
%     L = s;
%     elim = find(b+s>=0);
%     U = s;
%     U(elim) = b(elim);
%     eb(i) = errorbar(bins,b,L,U,symb{i},...
%         'color',cl(i,:),...
%         'linewidth',1.5,...
%         'markersize',8,...
%         'markeredgecolor','k',...
%         'markerfacecolor',cl(i,:));hold on
% end
% set(sp(1),'ylim',[0 0.03])
% set(sp(2),'ylim',[-0.03 0])
% title(sp(1),'Accretion')
% title(sp(2),'Erosion')
% xlabel(sp(2),'Event Length [min]')
% ylabel(sp(1),'Net bed change over event')
% ylabel(sp(2),'Net bed change over event')
% legend(eb,fn);
% %%%
% eb = zeros(3,1);
% for i = 1:length(fn)
%     %Try plotting uorb, umag, eventl against the bed movement stats
%     %Wave orbital velocity
%     f5 = figure(5);
%     set(f5,'PaperOrientation','portrait',...
%         'position',[400 100   800   400],...
%         'renderer','painters')
%     uorb = dat.(fn{i}).uorb;
%     deltbd = dat.(fn{i}).deltbd;
%     sp(1) = subplot(211);
%     id = find(deltbd>0);
%     bins = linspace(min(uorb(id)),max(uorb(id)),10);
%     [b,~,s] = bindata(uorb,deltbd(id),bins);
%     U = s;
%     elim = find(b-s<=0);
%     L = s;
%     L(elim) = b(elim);
%     eb(i) = errorbar(bins,b,L,U,symb{i},...
%         'color',cl(i,:),...
%         'linewidth',1.5,...
%         'markersize',8,...
%         'markeredgecolor','k',...
%         'markerfacecolor',cl(i,:));hold on
%     sp(2) = subplot(212);
%     id = find(deltbd<0);
%     bins = linspace(min(uorb(id)),max(uorb(id)),10);
%     [b,~,s] = bindata(uorb,deltbd(id),bins);
%     L = s;
%     elim = find(b+s>=0);
%     U = s;
%     U(elim) = b(elim);
%     eb(i) = errorbar(bins,b,L,U,symb{i},...
%         'color',cl(i,:),...
%         'linewidth',1.5,...
%         'markersize',8,...
%         'markeredgecolor','k',...
%         'markerfacecolor',cl(i,:));hold on
% end
% set(sp(1),'ylim',[0 0.03])
% set(sp(2),'ylim',[-0.03 0])
% title(sp(1),'Accretion')
% title(sp(2),'Erosion')
% xlabel(sp(2),'Wave Orbital Velocity [m/s]')
% ylabel(sp(1),'Net bed change over event')
% ylabel(sp(2),'Net bed change over event')
% legend(eb,fn);
% %%%
% eb = zeros(3,1);
% for i = 1:length(fn)
%     %Try plotting uorb, umag, eventl against the bed movement stats
%     %Cross shore velocity magnitude
%     f6 = figure(6);
%     set(f6,'PaperOrientation','portrait',...
%         'position',[400 100   800   400],...
%         'renderer','painters')
%     umag = dat.(fn{i}).umag;
%     deltbd = dat.(fn{i}).deltbd;
%     sp(1) = subplot(211);
%     id = find(deltbd>0);
%     bins = linspace(min(umag(id)),max(umag(id)),10);
%     [b,~,s] = bindata(umag,deltbd(id),bins);
%     U = s;
%     elim = find(b-s<=0);
%     L = s;
%     L(elim) = b(elim);
%     eb(i) = errorbar(bins,b,L,U,symb{i},...
%         'color',cl(i,:),...
%         'linewidth',1.5,...
%         'markersize',8,...
%         'markeredgecolor','k',...
%         'markerfacecolor',cl(i,:));hold on
%     sp(2) = subplot(212);
%     id = find(deltbd<0);
%     bins = linspace(min(umag(id)),max(umag(id)),10);
%     [b,~,s] = bindata(umag,deltbd(id),bins);
%     L = s;
%     elim = find(b+s>=0);
%     U = s;
%     U(elim) = b(elim);
%     eb(i) = errorbar(bins,b,L,U,symb{i},...
%         'color',cl(i,:),...
%         'linewidth',1.5,...
%         'markersize',8,...
%         'markeredgecolor','k',...
%         'markerfacecolor',cl(i,:));hold on
% end
% set(sp(1),'ylim',[0 0.03])
% set(sp(2),'ylim',[-0.03 0])
% title(sp(1),'Accretion')
% title(sp(2),'Erosion')
% xlabel(sp(2),'Cross-shore velocity magnitude [m/s]')
% ylabel(sp(1),'Net bed change over event')
% ylabel(sp(2),'Net bed change over event')
% legend(eb,fn);
%%%
%% Normalized Bed Variance
%%%
eb = zeros(3,1);
for i = 1:length(fn)
    %Try plotting uorb, umag, eventl against the bed movement stats
    %event length
    f7 = figure(7);
    set(f7,'PaperOrientation','portrait',...
        'position',[400 100   800   400],...
        'renderer','painters')
    eventl = dat.(fn{i}).eventl;
    bdvar = dat.(fn{i}).bdvar;
    bins = linspace(min(eventl),max(eventl),10);
    sp(1) = subplot(211);
    [b,~,s] = bindata(eventl,bdvar,bins);
    U = s;
    elim = find(b-s<=0);
    L = s;
    L(elim) = b(elim);
    pp = ploterr(bins,b,[],{L,U},symb{i},'hhx',0.01);
    eb(i) = pp(1);
    set(pp,...
        'color',cl(i,:),...
        'linewidth',1.5,...
        'markersize',8,...
        'markeredgecolor','k',...
        'markerfacecolor',cl(i,:));hold on
    sp(2) = subplot(212);
    [n1,x1] = hist(eventl,bins);
    bar(x1,n1./numel(eventl),...
        'facecolor',cl(i,:),...
        'edgecolor','none'),hold on
end
set(sp(1),'yscale','log')
xlabel('Event Length [min]')
ytext = '\sigma_{bd}/\Deltat*U [m]';
ylabel(sp(2),'% occurence')
ylabel(sp(1),ytext)
legend(eb,fn);
%%%
eb = zeros(3,1);
for i = 1:length(fn)
    %Try plotting uorb, umag, eventl against the bed movement stats
    %Wave Orbital Velocity
    f8 = figure(8);
    set(f8,'PaperOrientation','portrait',...
        'position',[400 100   800   400],...
        'renderer','painters')
    uorb = dat.(fn{i}).uorb;
    bdvar = dat.(fn{i}).bdvar;
    bins = linspace(min(uorb),max(uorb),10);
    sp(1) = subplot(211);
    [b,~,s] = bindata(uorb,bdvar,bins);
    U = s;
    elim = find(b-s<=0);
    L = s;
    L(elim) = b(elim);
    pp = ploterr(bins,b,[],{L,U},symb{i},'hhx',0.01);
    eb(i) = pp(1);
    set(pp,...
        'color',cl(i,:),...
        'linewidth',1.5,...
        'markersize',8,...
        'markeredgecolor','k',...
        'markerfacecolor',cl(i,:));hold on
    sp(2) = subplot(212);
    [n1,x1] = hist(uorb,bins);
    bar(x1,n1./numel(uorb),...
        'facecolor',cl(i,:),...
        'edgecolor','none'),hold on
end
set(sp(1),'yscale','log')
xlabel('Wave Orbital Velocity [m/s]')
ytext = '\sigma_{bd}/\Deltat*U [m]';
ylabel(sp(2),'% occurence')
ylabel(sp(1),ytext)
legend(eb,fn);
%%%
eb = zeros(3,1);
for i = 1:length(fn)
    %Try plotting uorb, umag, eventl against the bed movement stats
    %Cross shore velocity magnitude
    f9 = figure(9);
    set(f9,'PaperOrientation','portrait',...
        'position',[400 100   800   400],...
        'renderer','painters')
    umag = dat.(fn{i}).umag;
    bdvar = dat.(fn{i}).bdvar;
    bins = linspace(min(umag),max(umag),10);
    sp(1) = subplot(211);
    [b,~,s] = bindata(umag,bdvar,bins);
    U = s;
    elim = find(b-s<=0);
    L = s;
    L(elim) = b(elim);
    pp = ploterr(bins,b,[],{L,U},symb{i},'hhx',0.01);
    eb(i) = pp(1);
    set(pp,...
        'color',cl(i,:),...
        'linewidth',1.5,...
        'markersize',8,...
        'markeredgecolor','k',...
        'markerfacecolor',cl(i,:));hold on
    sp(2) = subplot(212);
    [n1,x1] = hist(umag,bins);
    bar(x1,n1./numel(umag),...
        'facecolor',cl(i,:),...
        'edgecolor','none'),hold on
end
set(sp(1),'yscale','log')
xlabel('Cross-shore velocity magnitude [m/s]')
ytext = '\sigma_{bd}/\Deltat*U [m]';
ylabel(sp(2),'% occurence')
ylabel(sp(1),ytext)
legend(eb,fn);
%%%
%Save figures
sfdir = 'f:\GradSchool\DataAnalysis\Paper3\WorkingFigures\Wavelet\';
prettyfigures('text',12,'labels',13,'box',1)
handles = findobj('type','figure');
label2 = {'umag';'uorb';'eventl'};label2 = repmat(label2,2,1);
for i = 1:3
    export_fig(figure(handles(i)),[sfdir 'CmbData_bdvar_' label2{i}],'-png')
end
for i = 4:6
    export_fig(figure(handles(i)),[sfdir 'CmbData_bdmed_' label2{i}],'-png')
end

