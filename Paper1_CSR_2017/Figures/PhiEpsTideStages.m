%Plot routine to generate cross-shore, phi/number of stems/stem diameter
%plots, plus phi,epsilon plots for various tide stages.

%make sure 'CatTurbData_v2.mat' is in the path

clear
%load the vegdat files and name based on the size and stage of tide
vpdir = 'd:\Projects\Mekong_W2015\DataAnalysis\DataReports\Paper1_QuadAnalysis\ThirdAttempt\';
vd = struct(); %big vegdat structure
order = [2 4 3 1]; %load files in order

%%%10cm%%%
files = dir([vpdir '*local10cm*.mat']);files = {files.name};
for i = 1:length(files)
    load([vpdir files{order(i)}])
    stage = regexp(files{order(i)},'.+_(.*).mat','tokens');
    fname = ['ten' char(stage{:})];
    vd.(fname).n = vegdat.n;
    vd.(fname).meanD = vegdat.meanD;
    vd.(fname).a = vegdat.a;
    vd.(fname).phi = vegdat.phi;
end
%%%20cm%%%
files = dir([vpdir '*local20cm*.mat']);files = {files.name};
for i = 1:length(files)
    load([vpdir files{order(i)}])
    stage = regexp(files{order(i)},'.+_(.*).mat','tokens');
    fname = ['twe' char(stage{:})];
    vd.(fname).n = vegdat.n;
    vd.(fname).meanD = vegdat.meanD;
    vd.(fname).a = vegdat.a;
    vd.(fname).phi = vegdat.phi;
end
%%%30cm%%%
files = dir([vpdir '*local30cm*.mat']);files = {files.name};
for i = 1:length(files)
    load([vpdir files{order(i)}])
    stage = regexp(files{order(i)},'.+_(.*).mat','tokens');
    fname = ['thr' char(stage{:})];
    vd.(fname).n = vegdat.n;
    vd.(fname).meanD = vegdat.meanD;
    vd.(fname).a = vegdat.a;
    vd.(fname).phi = vegdat.phi;
end
%%%40cm%%%
files = dir([vpdir '*local40cm*.mat']);files = {files.name};
for i = 1:length(files)
    load([vpdir files{order(i)}])
    stage = regexp(files{order(i)},'.+_(.*).mat','tokens');
    fname = ['fou' char(stage{:})];
    vd.(fname).n = vegdat.n;
    vd.(fname).meanD = vegdat.meanD;
    vd.(fname).a = vegdat.a;
    vd.(fname).phi = vegdat.phi;
end
%%%50cm%%%
files = dir([vpdir '*local50cm*.mat']);files = {files.name};
for i = 1:length(files)
    load([vpdir files{order(i)}])
    stage = regexp(files{order(i)},'.+_(.*).mat','tokens');
    fname = ['fif' char(stage{:})];
    vd.(fname).n = vegdat.n;
    vd.(fname).meanD = vegdat.meanD;
    vd.(fname).a = vegdat.a;
    vd.(fname).phi = vegdat.phi;
end
%%%1m%%%
load([vpdir 'Vegdat_1m.mat'])
vd.full.n = vegdat.n;
vd.full.meanD = vegdat.meanD;
vd.full.a = vegdat.a;
vd.full.phi = vegdat.phi;
vd.X = vegdat.Xshore;
clear vegdat
vegdat = vd;

%%%Velocity Statistics%%%
CatTurbData_v2
fn = fieldnames(vegdat);
E = [veldat.four.E; veldat.five.E];
std = [veldat.four.Estd; veldat.five.Estd];
H = [veldat.four.depth; veldat.five.depth];
Hs = [veldat.four.wrms; veldat.five.wrms].*sqrt(2);
gamma = Hs./H;
savefigdir = 'd:\Projects\Mekong_W2015\Figures\Paper2\Phi-Eps\';

%%%Plot Type%%%
lineplots = 0;
aggregate = 0;
csr = 1;
environment = 0;

if lineplots
    %%%Plot 10cm data%%%
    step = 0.005;
    f1 = figure(1);
    set(f1,'PaperOrientation','portrait',...
        'position',[400 100   800 600]);
    set(gcf,'color','w','PaperPositionMode','auto')
    c = brewermap(4,'Paired');
    
    colormap(c)
    vid = strfind(fn,'ten');
    id = find(not(cellfun('isempty',vid)));
    p = zeros(1,4);
    for i = 1:length(id)
        phi = vegdat.(fn{id(i)}).phi;
        zid = find(phi == 0);
        zid = setxor(1:54,zid);
        Ez = E(zid,i);Ez(Ez > 5E-3) = NaN;%Ez(Ez < 1E-4) = NaN;
        phiz = phi(zid);
        xs = min(phiz):step:max(phiz);
        [b,n,s] = bindata(phiz,Ez,xs);
        nanid = find(isnan(b));
        %remove NaNs
        xs = xs(setxor(1:length(xs),nanid));
        b = b(setxor(1:length(b),nanid));
        s = s(setxor(1:length(s),nanid));
        plot(xs,b,'-','LineWidth',1.5,...
            'Color',c(i,:)), hold on
        p(i) = errorbar(xs,b,s,'o',...
            'Color','k',...
            'LineWidth',1.5,...
            'markerfacecolor',c(i,:),...
            'markeredgecolor','k',...
            'markersize',10);
    end
    leg = legend(p,{'LL','ML','MH','HH'});
    set(leg,'box','off')
    axis([-0.025 1 -1E-4 1.2E-3])
    title('10cm^2 Sample Area','FontSize',14,'FontName','Cambria')
    xlabel('\phi','FontSize',14,'FontName','Cambria')
    ylabel('\epsilon (Wkg^-^1)','FontSize',14,'FontName','Cambria')
    set(gca,'LineWidth',1.5,'FontSize',13,'FontName','Cambria')
%     export_fig([savefigdir '10cm'],'-jpeg','-nocrop')
    
    %%%Plot 20cm data%%%
    f2 = figure(2);
    set(f2,'PaperOrientation','portrait',...
        'position',[400 100   800 600]);
    set(gcf,'color','w','PaperPositionMode','auto')
    c = brewermap(4,'Paired');
    
    colormap(c)
    vid = strfind(fn,'twe');
    id = find(not(cellfun('isempty',vid)));
    p = zeros(1,4);
    for i = 1:length(id)
        phi = vegdat.(fn{id(i)}).phi;
        zid = find(phi == 0);
        zid = setxor(1:54,zid);
        Ez = E(zid,i);Ez(Ez > 1E-3) = NaN;
        phiz = phi(zid);
        xs = min(phiz):step:max(phiz);
        [b,n,s] = bindata(phiz,Ez,xs);
        nanid = find(isnan(b));
        %remove NaNs
        xs = xs(setxor(1:length(xs),nanid));
        b = b(setxor(1:length(b),nanid));
        s = s(setxor(1:length(s),nanid));
        plot(xs,b,'-','LineWidth',1.5,...
            'Color',c(i,:)), hold on
        p(i) = errorbar(xs,b,s,'o',...
            'Color','k',...
            'LineWidth',1.5,...
            'markerfacecolor',c(i,:),...
            'markeredgecolor','k',...
            'markersize',10);
    end
    leg = legend(p,{'LL','ML','MH','HH'});
    set(leg,'box','off')
    axis([-0.025 1 -1E-4 1.2E-3])
    title('20cm^2 Sample Area','FontSize',14,'FontName','Cambria')
    xlabel('\phi','FontSize',14,'FontName','Cambria')
    ylabel('\epsilon (Wkg^-^1)','FontSize',14,'FontName','Cambria')
    set(gca,'LineWidth',1.5,'FontSize',13,'FontName','Cambria')
%     export_fig([savefigdir '20cm'],'-jpeg','-nocrop')
    
    %%%Plot 30cm data%%%
    step = 0.005;
    f3 = figure(3);
    set(f3,'PaperOrientation','portrait',...
        'position',[400 100   800 600]);
    set(gcf,'color','w','PaperPositionMode','auto')
    c = brewermap(4,'Paired');
    
    colormap(c)
    vid = strfind(fn,'thr');
    id = find(not(cellfun('isempty',vid)));
    p = zeros(1,4);
    for i = 1:length(id)
        phi = vegdat.(fn{id(i)}).phi;
        zid = find(phi == 0);
        zid = setxor(1:54,zid);
        Ez = E(zid,i);Ez(Ez > 1E-3) = NaN;
        phiz = phi(zid);
        xs = min(phiz):step:max(phiz);
        [b,n,s] = bindata(phiz,Ez,xs);
        nanid = find(isnan(b));
        %remove NaNs
        xs = xs(setxor(1:length(xs),nanid));
        b = b(setxor(1:length(b),nanid));
        s = s(setxor(1:length(s),nanid));
        plot(xs,b,'-','LineWidth',1.5,...
            'Color',c(i,:)), hold on
        p(i) = errorbar(xs,b,s,'o',...
            'Color','k',...
            'LineWidth',1.5,...
            'markerfacecolor',c(i,:),...
            'markeredgecolor','k',...
            'markersize',10);
    end
    leg = legend(p,{'LL','ML','MH','HH'});
    set(leg,'box','off')
    axis([-0.025 0.5 -1E-4 1.2E-3])
    title('30cm^2 Sample Area','FontSize',14,'FontName','Cambria')
    xlabel('\phi','FontSize',14,'FontName','Cambria')
    ylabel('\epsilon (Wkg^-^1)','FontSize',14,'FontName','Cambria')
    set(gca,'LineWidth',1.5,'FontSize',13,'FontName','Cambria')
%     export_fig([savefigdir '30cm'],'-jpeg','-nocrop')
    
    %%%Plot 40cm data%%%
    step = 0.005;
    f4 = figure(4);
    set(f4,'PaperOrientation','portrait',...
        'position',[400 100   800 600]);
    set(gcf,'color','w','PaperPositionMode','auto')
    c = brewermap(4,'Paired');
    
    colormap(c)
    vid = strfind(fn,'fou');
    id = find(not(cellfun('isempty',vid)));
    p = zeros(1,4);
    for i = 1:length(id)
        phi = vegdat.(fn{id(i)}).phi;
        zid = find(phi == 0);
        zid = setxor(1:54,zid);
        Ez = E(zid,i);Ez(Ez > 5E-3) = NaN;Ez(Ez < 1E-4) = NaN;
        phiz = phi(zid);
        xs = min(phiz):step:max(phiz);
        [b,n,s] = bindata(phiz,Ez,xs);
        nanid = find(isnan(b));
        %remove NaNs
        xs = xs(setxor(1:length(xs),nanid));
        b = b(setxor(1:length(b),nanid));
        s = s(setxor(1:length(s),nanid));
        plot(xs,b,'-','LineWidth',1.5,...
            'Color',c(i,:)), hold on
        p(i) = errorbar(xs,b,s,'o',...
            'Color','k',...
            'LineWidth',1.5,...
            'markerfacecolor',c(i,:),...
            'markeredgecolor','k',...
            'markersize',10);
    end
    leg = legend(p,{'LL','ML','MH','HH'});
    set(leg,'box','off')
    axis([-0.025 0.5 -1E-4 1.2E-3])
    title('40cm^2 Sample Area','FontSize',14,'FontName','Cambria')
    xlabel('\phi','FontSize',14,'FontName','Cambria')
    ylabel('\epsilon (Wkg^-^1)','FontSize',14,'FontName','Cambria')
    set(gca,'LineWidth',1.5,'FontSize',13,'FontName','Cambria')
%     export_fig([savefigdir '40cm'],'-jpeg','-nocrop')
    
    %%%Plot 50cm data%%%
    f5 = figure(5);
    set(f5,'PaperOrientation','portrait',...
        'position',[400 100   800 600]);
    set(gcf,'color','w','PaperPositionMode','auto')
    c = brewermap(4,'Paired');
    
    colormap(c)
    vid = strfind(fn,'fif');
    id = find(not(cellfun('isempty',vid)));
    p = zeros(1,4);
    for i = 1:length(id)
        phi = vegdat.(fn{id(i)}).phi;
        zid = find(phi == 0);
        zid = setxor(1:54,zid);
        Ez = E(zid,i);Ez(Ez > 5E-3) = NaN;
        phiz = phi(zid);
        xs = min(phiz):step:max(phiz);
        [b,n,s] = bindata(phiz,Ez,xs);
        nanid = find(isnan(b));
        %remove NaNs
        xs = xs(setxor(1:length(xs),nanid));
        b = b(setxor(1:length(b),nanid));
        s = s(setxor(1:length(s),nanid));
        plot(xs,b,'-','LineWidth',1.5,...
            'Color',c(i,:)), hold on
        p(i) = errorbar(xs,b,s,'o',...
            'Color','k',...
            'LineWidth',1.5,...
            'markerfacecolor',c(i,:),...
            'markeredgecolor','k',...
            'markersize',10);
    end
    leg = legend(p,{'LL','ML','MH','HH'});
    set(leg,'box','off')
    axis([-0.025 0.5 -1E-4 1.2E-3])
    title('50cm^2 Sample Area','FontSize',14,'FontName','Cambria')
    xlabel('\phi','FontSize',14,'FontName','Cambria')
    ylabel('\epsilon (Wkg^-^1)','FontSize',14,'FontName','Cambria')
    set(gca,'LineWidth',1.5,'FontSize',13,'FontName','Cambria')
%     export_fig([savefigdir '50cm'],'-jpeg','-nocrop')
    
    %%%Plot 1m data%%%
    step = 5E-6;
    f6 = figure(6);
    set(f6,'PaperOrientation','portrait',...
        'position',[400 100   800 600]);
    set(gcf,'color','w','PaperPositionMode','auto')
    phi = vegdat.full.phi;
    Ez = mean(E,2);Ez(Ez > 1E-3) = NaN;
    phiz = phi;
    xs = min(phiz):step:max(phiz);
    [b,~,s] = bindata(phiz,Ez,xs);
    nanid = find(isnan(b));
    %remove NaNs
    xs = xs(setxor(1:length(xs),nanid));
    b = b(setxor(1:length(b),nanid));
    s = s(setxor(1:length(s),nanid));
    plot(xs,b,'-','LineWidth',1.5,...
        'Color','r'), hold on
    errorbar(xs,b,s,'o',...
        'Color','k',...
        'LineWidth',1.5,...
        'markerfacecolor','r',...
        'markeredgecolor','k',...
        'markersize',10);
    axis([-0.001 6E-3 -1E-4 1E-3])
    set(gca,'LineWidth',1.5,...
        'FontSize',14,...
        'FontName','Arial',...
        'TickDir','out')
    title('1m^2 Sample Area','FontSize',14,'FontName','Arial')
    xlabel('\phi','FontSize',18,'FontName','Arial')
    ylabel('\epsilon (Wkg^-^1)','FontSize',18,'FontName','Arial')
%     export_fig([savefigdir '1m'],'-jpeg','-nocrop')
    
    f7 = figure(7);
    set(f7,'PaperOrientation','portrait',...
        'position',[400 100   800 600]);
    set(gcf,'color','w','PaperPositionMode','auto')
    phi = vegdat.full.phi;
    Ez = mean(E,2);Ez(Ez > 1E-3) = NaN;
    loglog(phi,Ez,'o',...
        'Color','k',...
        'LineWidth',1.5,...
        'MarkerFaceColor','r',...
        'MarkerSize',10);
    set(gca,'LineWidth',1.5,...
        'FontSize',14,...
        'FontName','Arial',...
        'TickDir','out')
    title('1m^2 Sample Area - All Samples','FontSize',14,'FontName','Arial')
    xlabel('\phi','FontSize',18,'FontName','Arial')
    ylabel('\epsilon (Wkg^-^1)','FontSize',18,'FontName','Arial')
%     export_fig([savefigdir '1m_allsamples'],'-jpeg','-nocrop')

end
if aggregate
%     %%%Plot the whole dataset%%%
%     symb = {'^','o','d','s','<'};
%     tt = {'LL','ML','MH','HH'};
%     cb = {'YlOrRd';'PuBu';'BuPu';'YlGn'};
%     f1 = figure(1);
%     set(f1,'PaperOrientation','portrait',...
%         'position',[400 100   1000 800]);
%     set(gcf,'color','w','PaperPositionMode','auto')
%     sp = zeros(1,4);
%     p = zeros(1,5);
%     for i = 1:4
%         c = brewermap(5,cb{i});
% %         sp(i) = subplot(2,2,i);
%         vid = strfind(fn,tt{i});
%         id = find(not(cellfun('isempty',vid)));
%         for j = 1:length(id)
%             phi = vegdat.(fn{id(j)}).phi;
%             a = vegdat.(fn{id(j)}).a;
%             zid = find(a == 0);
%             zid = setxor(1:54,zid);
%             Ez = E(zid,i);
%             az = a(zid);
% %             xs = min(phiz):step(i,j):max(phiz);
% %             [b,n,s] = bindata(phiz,Ez,xs);
% %             nanid = find(isnan(b));
% %             %remove NaNs
% %             xs = xs(setxor(1:length(xs),nanid));
% %             b = b(setxor(1:length(b),nanid));
% %             s = s(setxor(1:length(s),nanid));
%             plot(az,Ez,symb{j},...
%                 'markersize',10,...
%                 'markerfacecolor',c(j,:),...
%                 'markeredgecolor','k'), hold on
% %             plot(xs,b,'-','LineWidth',1.5,...
% %                 'Color',c(j,:)), hold on
% %             p(j) = errorbar(xs,b,s,symb{j},...
% %                 'Color','k',...
% %                 'LineWidth',1.5,...
% %                 'markerfacecolor',c(j,:),...
% %                 'markeredgecolor','k',...
% %                 'markersize',10);
%         end
%         text(0.85,9E-4,tt{i},'FontSize',18,...
%             'FontName','Arial')
%     end
%     set(gca,'LineWidth',1.5,...
%         'FontSize',14,...
%         'FontName','Arial',...
%         'Xlim',[-0.1 9],...
%         'YScale','log',...
%         'TickDir','out')
%     set(sp(1),'position',[0.1 0.54 0.4 0.4],...
%         'XTickLabel',[])
%     set(sp(2),'position',[0.53 0.54 0.4 0.4],...
%         'YTickLabel',[],'XTickLabel',[])
%     set(sp(3),'position',[0.1 0.1 0.4 0.4])
%     set(sp(4),'position',[0.53 0.1 0.4 0.4],...
%         'YTickLabel',[])
%     xlabel(sp(3),'\phi','FontSize',18,'FontName','Arial')
%     xlabel(sp(4),'\phi','FontSize',18,'FontName','Arial')
%     ylabel(sp(1),'\epsilon (Wkg^-^1)','FontSize',18,'FontName','Arial')
%     ylabel(sp(3),'\epsilon (Wkg^-^1)','FontSize',18,'FontName','Arial')
%     leg = legend(p,{'10cm^2';'20cm^2';'30cm^2';'40cm^2';'50cm^2'});
%     set(leg,'position',[0.85 0.32 0.08 0.08],...
%         'LineWidth',1.5)
%         export_fig([savefigdir 'AllWindows_vol'],'-pdf','-nocrop')
% 
%     %%Plot an example 'scattershot' of Phi/Eps
%     tt = {'ten','twe','thr','fou','fif','full'};
%     f2 = figure(2);
%     set(f2,'PaperOrientation','portrait',...
%         'position',[400 100   600 500]);
%     set(gcf,'color','w','PaperPositionMode','auto')
%     p = zeros(1,5);
%     step = 0.001;
%     phi = zeros(54,21);
%     ct = 1;
%     for i = 1:6
%         vid = strfind(fn,tt{i});
%         id = find(not(cellfun('isempty',vid)));
%         for j = 1:length(id)
%             phi(:,ct) = vegdat.(fn{id(j)}).phi;
%             ct = ct+1;
%         end
%     end
%     phi = mean(phi,2);zid = find(phi == 0);
%     zid = setxor(1:54,zid);
%     Ez = mean(E(zid,:),2);
%     phiz = phi(zid);
%     xs = min(phiz):step:max(phiz);
%     [b,n,s] = bindata(phiz,Ez,xs);
%     errorbar(xs,log10(b),log10(s),'o',...
%         'Color','k',...
%         'LineWidth',1.5,...
%         'markerfacecolor',[0.8 0.8 0.8],...
%         'markeredgecolor','k',...
%         'markersize',10);
%     xlabel('\phi','FontSize',18,'FontName','Arial')
%     ylabel('log_1_0(\epsilon) (Wkg^-^1)','FontSize',18,'FontName','Arial')
%     set(gca,...
%         'XLim',[-0.001 0.025],...
%         'YLim',[-9 2],...
%         'YTick',-8:2:2,...
%         'FontSize',14,...
%         'FontName','Arial',...
%         'LineWidth',1.5,...
%         'TickDir','out')
%     export_fig([savefigdir 'Xphi_noscale'],'-pdf','-nocrop')
    
end
if csr
    %%%X and a for a single scale (20cm^2)
%     step = 10;
%     f1 = figure(1);
%     set(f1,'PaperOrientation','portrait',...
%         'position',[400 100   1000 500]);
%     set(gcf,'color','w','PaperPositionMode','auto')
%     vid = strfind(fn,'twe');
%     id = find(not(cellfun('isempty',vid)));
%     
%     sp(1) = subplot(211);
%     a = zeros(54,4);
%     for i = 1:length(id)
%         a(:,i) = vegdat.(fn{id(i)}).a;
%     end
%     a = nanmean(a,2);
%     X = vegdat.X;
%     xs = min(X):step:max(X);
%     [b,~,s] = bindata(X,a,xs);
%     nanid = find(isnan(b));
%     %remove NaNs
%     xs = xs(setxor(1:length(xs),nanid));
%     b = b(setxor(1:length(b),nanid));
%     s = s(setxor(1:length(s),nanid));
% %     s(5)= 0.016;
%     plot(xs,b,'-','LineWidth',1.5,...
%         'Color','k'), hold on
%     errorbar(xs,b,s,'o',...
%         'Color','k',...
%         'LineWidth',1.5,...
%         'markerfacecolor',[0.8 0.8 0.8],...
%         'markeredgecolor','k',...
%         'markersize',10);
%     set(gca,'LineWidth',1.5,...
%         'FontSize',13,...
%         'FontName','Cambria',...
%         'YLim',[-0.005 0.15],...
%         'YTick',0:0.03:0.15,...
%         'TickDir','out')
%     ylabel('a (m^-^1)','FontSize',14)
%     
%     sp(2) = subplot(212);
%     Ez = nanmean(E,2);Ez(Ez > 1E-3) = NaN;
%     [b,~,s] = bindata(X,Ez,xs);
%     nanid = find(isnan(b));
%     %remove NaNs
%     xs = xs(setxor(1:length(xs),nanid));
%     b = b(setxor(1:length(b),nanid));
%     s = s(setxor(1:length(s),nanid));
%     plot(xs,b,'-','LineWidth',1.5,...
%         'Color','k'), hold on
%     errorbar(xs,b,s,'^',...
%         'Color','k',...
%         'LineWidth',1.5,...
%         'markerfacecolor',[0.8 0.8 0.8],...
%         'markeredgecolor','k',...
%         'markersize',10);
%     set(gca,'LineWidth',1.5,...
%         'FontSize',13,...
%         'FontName','Cambria',...
%         'YLim',[-1E-4 2E-3],...
%         'TickDir','out')
%     ylabel('\epsilon (Wkg^-^1)','FontSize',14)
%     xlabel('Cross-shore distance (m)','FontSize',14)
%     
%     set(sp(1),'Position',[0.1 0.55 0.84 0.37],...
%         'XTickLabel',[])
%     set(sp(2),'Position',[0.1 0.12 0.84 0.37])
%     hold off
%         export_fig([savefigdir 'Xa_Xeps'],'-jpeg','-nocrop')

    %Plot wave breaking (gamma) and epsilon for all four tidal stages
    symb = {'^','o','d','s'};
    f1 = figure(1);
    set(f1,'PaperOrientation','portrait',...
        'position',[400 100   800 500]);
    set(gcf,'color','w','PaperPositionMode','auto')
    c = brewermap(4,'YlGnBu');colormap(c)
    c(1,:) = [255 235 0]./255;
    p = zeros(1,4);
    for i = 1:4
        p(i) = plot(gamma(:,i),E(:,i),symb{i},...
            'LineWidth',1.5,...
            'MarkerSize',9,...
            'MarkerFaceColor',c(i,:),...
            'MarkerEdgeColor','k'); hold on
    end
   leg = legend(p,{'LL';'ML';'MH';'HH'});
   set(leg,'position',[0.87 0.26 0.005 0.005])
   set(gca,...
       'Yscale','log',...
       'Box','on',...
       'LineWidth',1.5,...
       'FontSize',14,...
       'FontName','Arial',...
       'TickDir','out',...
       'XLim',[0 0.6])
   ylabel('\epsilon (Wkg^-^1)','Fontsize',18)
   xlabel('\gamma','FontSize',18)
%    export_fig([savefigdir 'GammaEps'],'-pdf','-nocrop')

    %%%X and a for a single scale, including LL to HH (20cm^2)
    symb = {'^','o','d','s'};
    lins = {':','-.','--','-'};
    f2 = figure(2);
    set(f2,'PaperOrientation','portrait',...
        'position',[400 100   1000 800]);
    set(gcf,'color','w','PaperPositionMode','auto')
    sp(1) = subplot(311);

    N = zeros(54,5);p = zeros(1,5);
    Hs = nanmean(Hs,2);
    xs = -60:10:100;
    [b,n,s] = bindata(vegdat.X,Hs,xs);
    nanid = find(isnan(b));
    %remove NaNs
    xs = xs(setxor(1:length(xs),nanid));
    b = b(setxor(1:length(b),nanid));
    s = s(setxor(1:length(s),nanid));
    line([0 0],[0 1],...
        'LineStyle','--',...
        'LineWidth',1.5,...
        'Color','k'), hold on
    errorbar(xs,b,s,'d','Color','k',...
        'LineWidth',1.5,...
        'markerfacecolor',[0,206,209]./255,...
        'markeredgecolor','k',...
        'markersize',12)
    axis([-80 100 0 0.5])
    set(gca,'LineWidth',1.5,'box','on',...
        'FontSize',14,'FontName','Arial',...
        'TickDir','out',...
        'TickLength',[0.009 0.009],...
        'YTick',0:0.1:0.5,...
        'XTickLabel',[])
    ylabel('H_s (m)','FontSize',18)

    sp(2) = subplot(312);
    vid = strfind(fn,'twe');
    id = find(not(cellfun('isempty',vid)));
    step = 10;a = zeros(54,4);
    for i = 1:length(id)
        a(:,i) = vegdat.(fn{id(i)}).a;
    end
    a = nanmean(a,2);
    X = vegdat.X;
    xs = min(X):step:max(X);
    [b,~,s] = bindata(X,a,xs);
    nanid = find(isnan(b));
    %remove NaNs
    xs = xs(setxor(1:length(xs),nanid));
    b = b(setxor(1:length(b),nanid));
    s = s(setxor(1:length(s),nanid));
%     s(5)= 0.016;
    plot(xs,b,'-','LineWidth',2,...
        'Color','k'), hold on
    errorbar(xs,b,s,'o',...
        'Color','k',...
        'LineWidth',1.5,...
        'markerfacecolor',[0.8 0.8 0.8],...
        'markeredgecolor','k',...
        'markersize',10);
    set(gca,'LineWidth',1.5,...
        'FontSize',14,...
        'FontName','Arial',...
        'YLim',[-0.1 4],...
        'YTick',0:1:4,...
        'TickDir','out')
    ylabel('a (m^-^1)','FontSize',18)
    
    sp(3) = subplot(313);
    c = brewermap(4,'YlGnBu');colormap(c)
    c(1,:) = [255 235 0]./255;
    p = zeros(1,4);
    for i = 1:4
        Ez = E(:,i);
        [b,~,s] = bindata(X,Ez,xs);
        nanid = find(isnan(b));
        %remove NaNs
        xs = xs(setxor(1:length(xs),nanid));
        b = b(setxor(1:length(b),nanid));
        s = s(setxor(1:length(s),nanid));
        plot(xs,b,lins{i},'LineWidth',2,...
            'Color',c(i,:)), hold on
        p(i) = errorbar(xs,b,s,symb{i},...
            'Color','k',...
            'LineWidth',1.5,...
            'markerfacecolor',c(i,:),...
            'markeredgecolor','k',...
            'markersize',10);
%         plot(X,Ez,'o','color',c(i,:)), hold on
    end
    set(gca,'LineWidth',1.5,...
        'FontSize',14,...
        'FontName','Arial',...
        'TickDir','out',...
        'YScale','log')
    leg = legend(p,{'LL';'ML';'MH';'HH'});
    set(leg,'position',[0.92 0.24 0.05 0.05],...
        'LineWidth',1.5)
    ylabel('\epsilon (Wkg^-^1)','FontSize',18)
    xlabel('Cross-shore distance (m)','FontSize',18)
    
    set(sp(1),'Position',[0.1 0.68 0.84 0.23],...
        'XTickLabel',[])
    set(sp(2),'Position',[0.1 0.4 0.84 0.23],...
        'XTickLabel',[])
    set(sp(3),'Position',[0.1 0.12 0.84 0.23])
    hold off
% 	export_fig([savefigdir 'Xa_SP_4'],'-pdf','-nocrop')
    
    
    %plot a and eps, phi and eps. Find values of each that correspond to
    %breaking gamma > 0.3 and leave as 'unfilled' symbols
    
    symb = {'^','o','d','s'};
    lins = {':','-.','--','-'};
    f3 = figure(3);
    set(f3,'PaperOrientation','portrait',...
        'position',[400 100   1000 500]);
    set(gcf,'color','w','PaperPositionMode','auto')
    sp(1) = subplot(121);
    c = brewermap(4,'YlGnBu');colormap(c)
    c(1,:) = [255 235 0]./255;
    vid = strfind(fn,'twe');
    id = find(not(cellfun('isempty',vid)));
    p = zeros(1,4);
    for i = 1:length(id)
        a = vegdat.(fn{id(i)}).a;
        zid = find(a == 0);
        zid = setxor(1:54,zid);
        Ez = E(zid);
        az = a(zid);        
        gam = gamma(zid,i);
        gid = find(gam >= 0.3);
        gdx = setxor(1:length(az),gid);
        p(i) = plot(az(gdx),Ez(gdx),symb{i},...
            'LineWidth',1.5,...
            'Markersize',9,...
            'MarkerFaceColor',c(i,:),...
            'MarkerEdgeColor','k');hold on
        p(i) = plot(az(gid),Ez(gid),symb{i},...
            'LineWidth',1.5,...
            'Markersize',9,...
            'MarkerFaceColor','none',...
            'MarkerEdgeColor',c(i,:));        
    end
    set(gca,'LineWidth',1.5,...
        'YScale','log',...
        'FontSize',14,...
        'FontName','Arial',...
        'YLim',[-1E-4 1.2E-2],...
        'XLim',[-1E-1 7],...
        'XTick',0:1:7,...
        'TickDir','out')
    ylabel('\epsilon (Wkg^-^1)','FontSize',18)
    xlabel('a (m^-^1)','FontSize',18)
    
    sp(2) = subplot(122);
    p = zeros(1,4);
    for i = 1:length(id)
        phi = vegdat.(fn{id(i)}).phi;
        zid = find(phi == 0);
        zid = setxor(1:54,zid);
        Ez = E(zid);
        phiz = phi(zid);
        gam = gamma(zid,i);
        gid = find(gam >= 0.3);
        gdx = setxor(1:length(az),gid);
        p(i) = plot(phiz(gdx),Ez(gdx),symb{i},...
            'LineWidth',1.5,...
            'Markersize',9,...
            'MarkerFaceColor',c(i,:),...
            'MarkerEdgeColor','k');     
        hold on
        q(i) = plot(phiz(gid),Ez(gid),symb{i},...
            'LineWidth',1.5,...
            'Markersize',9,...
            'MarkerFaceColor','none',...
            'MarkerEdgeColor',c(i,:)); 
    end
    set(gca,'LineWidth',1.5,...
        'YScale','log',...
        'FontSize',14,...
        'FontName','Arial',...
        'YLim',[-1E-4 1.2E-2],...
        'YTickLabel',[],...
        'XLim',[-5E-3 0.1],...
        'XTick',0:0.02:0.1,...
        'TickDir','out')
    leg = legend(p,{'LL';'ML';'MH';'HH'});
    set(leg,'position',[0.9 0.26 0.05 0.05],...
        'LineWidth',1.5,'box','on')
    xlabel('\phi','FontSize',18)
    set(sp(1),'position',[0.11 0.15 0.4 0.8])
    set(sp(2),'position',[0.55 0.15 0.4 0.8])
    % 	export_fig([savefigdir 'A_Phi_eps_&Gamma_v2'],'-pdf','-nocrop')
    
    f4 = figure(4);
    set(f4,'PaperOrientation','portrait',...
        'position',[400 100   800 500]);
    set(gcf,'color','w','PaperPositionMode','auto')
    c = [0.2 0.2 0.2;0.4 0.4 0.4;0.6 0.6 0.6;0.8 0.8 0.8];
    subplot(121)
    vid = strfind(fn,'twe');
    id = find(not(cellfun('isempty',vid)));
    A = zeros(54,4);
    for i = 1
        A = vegdat.(fn{id(i)}).a;
        eps = E(:,i);
        A(isnan(eps)) = 0;
        eps(isnan(eps)) = 0;        
        zid = find(A == 0);
        zid = setxor(1:54,zid);
        eps = eps(zid);
        Az = A(zid);
        pf = polyfit(log(Az),log(eps),1);
        y_hat = exp(pf(1)*log(Az)+pf(2));
        loglog(Az,y_hat,'-',...
            'color',c(i,:),...
            'linewidth',1.5)
        hold on
        loglog(Az,eps,'+',...
            'color',c(i,:))
        disp(fn{id(i)})
        fitlm(log(Az),log(eps)) 
    end
    title('20 cm^2 Quadrats')
    xlabel('a (m^-^1)')
    ylabel('\epsilon (Wkg^-^1)')
    subplot(122)
    p = zeros(1,4);
    for i = 1
        A = vegdat.full.a;
        eps = E(:,i);
        A(isnan(eps)) = 0;
        eps(isnan(eps)) = 0;
        zid = find(A == 0);
        zid = setxor(1:54,zid);
        Az = A(zid);
        eps = eps(zid);
        pf = polyfit(log(Az),log(eps),1);
        y_hat = exp(pf(1)*log(Az)+pf(2));
        p(i) = loglog(Az,y_hat,'-',...
            'color',c(i,:),...
            'linewidth',1.5);
        hold on
        loglog(Az,eps,'*',...
            'color',c(i,:))
        disp(fn{id(i)})
        fitlm(log(Az),log(eps))
    end
    title('1 m^2 Quadrats')
    xlabel('a (m^-1)')
    leg = legend(p,{'LL';'ML';'MH';'HH'});
    prettyfigures('text',14,'labels',16,...
        'box',1)
    set(leg,'box','off')
    % export_fig([savefigdir '20cm_1m_compare'],'-jpeg','-nocrop')
    
    %plot x and eps categorized by the breaking condition
    symb = {'^','o','d','s'};
    f6 = figure(6);
    set(f6,'PaperOrientation','portrait',...
        'position',[400 100   1000 300]);
    set(gcf,'color','w','PaperPositionMode','auto')
    c = brewermap(4,'YlGnBu');colormap(c)
    c(1,:) = [255 235 0]./255;
    xs = -60:10:100;
    vid = strfind(fn,'twe');
    id = find(not(cellfun('isempty',vid)));
    p = zeros(1,4);
    for i = 1:length(id)
        Ez = E(:,i);
        gam = gamma(:,i);
        gid = find(gam >= 0.3);
        gdx = setxor(1:length(Ez),gid);
        %         p(i) = plot(X(gdx),Ez(gdx),symb{i},...
        %             'LineWidth',1.5,...
        %             'Markersize',9,...
        %             'MarkerFaceColor',c(i,:),...
        %             'MarkerEdgeColor','k');
        
        hold on
        %         p(i) = plot(X(gid),Ez(gid),symb{i},...
        %             'LineWidth',1.5,...
        %             'Markersize',9,...
        %             'MarkerFaceColor','none',...
        %             'MarkerEdgeColor',c(i,:));
        [b,~,s] = bindata(X(gid),Ez(gid),xs);
        nanid = find(isnan(b));
        %remove NaNs
        xs = xs(setxor(1:length(xs),nanid));
        b = b(setxor(1:length(b),nanid));
        s = s(setxor(1:length(s),nanid));
%         plot(xs,b,lins{i},'LineWidth',2,...
%             'Color',c(i,:)), hold on
        p(i) = errorbar(xs,b,s,symb{i},...
            'Color','k',...
            'LineWidth',1.5,...
            'markerfacecolor','w',...
            'markeredgecolor',c(i,:),...
            'markersize',10);
    end
    ylabel('\epsilon (Wkg^-^1)')
    xlabel('Cross-shore distance (m)')
    leg = legend(p,{'LL';'ML';'MH';'HH'});
    prettyfigures('text',14,'labels',16,...
        'box',1)
    set(leg,'position',[0.91 0.71 0.02 0.02])
    set(gca,'yscale','log')
    export_fig([savefigdir 'Xa_gamma'],'-pdf','-nocrop')
    
end
if environment
    %%%X,n%%%
    f1 = figure(1);
    set(f1,'PaperOrientation','portrait',...
        'position',[400 100   1200 700]);
    set(gcf,'color','w','PaperPositionMode','auto')
    tt = {'ten','twe','thr','fou','fif'};
    c = brewermap(5,'YlOrRd');colormap(c)
    xs = -60:5:100;
    [b,~,s] = bindata(vegdat.X,vegdat.full.n,xs);
    nanid = find(isnan(b));
    %remove NaNs
    xs = xs(setxor(1:length(xs),nanid));
    b = b(setxor(1:length(b),nanid));
    s = s(setxor(1:length(s),nanid));
    sp(1) = subplot(211);
    line([0 0],[0 200],...
        'LineStyle','--',...
        'LineWidth',1.5,...
        'Color','k'), hold on
    plot(xs,b,'LineWidth',1.5,'Color','k'), hold on
    p = errorbar(xs,b,s,'o','Color','k',...
        'LineWidth',1.5,...
        'markerfacecolor',[0.8 0.8 0.8],...
        'markeredgecolor','k',...
        'markersize',10);
    axis([-80 100 -10 140])
    set(gca,'LineWidth',1.5,'box','on',...
        'FontSize',13,'FontName','Cambria',...
        'XTickLabel',[],'TickDir','out',...
        'TickLength',[0.009 0.009])
    leg = legend(p,'1m');set(leg,'box','off')
    ylabel('Number of stems','FontSize',14)
    
    sp(2) = subplot(212);
    line([0 0],[0 200],...
        'LineStyle','--',...
        'LineWidth',1.5,...
        'Color','k'), hold on
    for i = 1:5
        vid = strfind(fn,tt{i});
        id = find(not(cellfun('isempty',vid)));
        n = zeros(54,4);
        for j = 1:length(id)
            n(:,j) = vegdat.(fn{id(j)}).n;
        end
        N = ceil(nanmean(n,2));
        xs = -60:10:100;
        [b,~,s] = bindata(vegdat.X,N,xs);
        nanid = find(isnan(b));
        %remove NaNs
        xs = xs(setxor(1:length(xs),nanid));
        b = b(setxor(1:length(b),nanid));
        s = s(setxor(1:length(s),nanid));
        plot(xs,b,'LineWidth',1.5,'Color',c(i,:)), hold on
        p(i) = errorbar(xs,b,s,'o',...
            'Color','k',...
            'LineWidth',1.5,...
            'markerfacecolor',c(i,:),...
            'markeredgecolor','k',...
            'markersize',10);
    end
    axis([-80 100 -2 40])
    set(gca,'LineWidth',1.5,'box','on',...
        'FontSize',13,'FontName','Cambria',...
        'TickDir','out','TickLength',[0.009 0.009])
    leg = legend(p,{'10cm^2';'20cm^2';'30cm^2';'40cm^2';'50cm^2'});
    set(leg,'box','off')
    ylabel('Number of stems','FontSize',14)
    xlabel('Cross-shore distance (m)','FontSize',14)
    set(sp(1),'position',[0.1 0.53 0.85 0.36])
    set(sp(2),'position',[0.1 0.1 0.85 0.36])
%     export_fig([savefigdir 'Xn'],'-jpeg','-nocrop')
    
    %%%X,gamma%%%
    f2 = figure(2);
    set(f2,'PaperOrientation','portrait',...
        'position',[400 100   1000 400]);
    set(gcf,'color','w','PaperPositionMode','auto')
    N = zeros(54,5);p = zeros(1,5);
    gam = nanmean(gamma,2);
    xs = -60:5:100;
    [b,n,s] = bindata(vegdat.X,gam,xs);
    nanid = find(isnan(b));
    %remove NaNs
    xs = xs(setxor(1:length(xs),nanid));
    b = b(setxor(1:length(b),nanid));
    s = s(setxor(1:length(s),nanid));
    line([0 0],[0 1],...
        'LineStyle','--',...
        'LineWidth',1.5,...
        'Color','k'), hold on
    errorbar(xs,b,s,'d','Color','k',...
        'LineWidth',1.5,...
        'markerfacecolor',[0,206,209]./255,...
        'markeredgecolor','k',...
        'markersize',12)
    axis([-80 100 0 0.7])
    set(gca,'LineWidth',1.5,'box','on',...
        'FontSize',14,'FontName','Arial',...
        'TickDir','out',...
        'TickLength',[0.009 0.009])
    ylabel('\gamma','FontSize',18)
    xlabel('Cross-shore distance (m)','FontSize',18)
%     export_fig([savefigdir 'Xgamma'],'-pdf','-nocrop')

    %%%X,Hs%%%
    f3 = figure(3);
    set(f3,'PaperOrientation','portrait',...
        'position',[400 100   1000 250]);
    set(gcf,'color','w','PaperPositionMode','auto')
    N = zeros(54,5);p = zeros(1,5);
    Hs = nanmean(Hs,2);
    xs = -60:10:100;
    [b,n,s] = bindata(vegdat.X,Hs,xs);
    nanid = find(isnan(b));
    %remove NaNs
    xs = xs(setxor(1:length(xs),nanid));
    b = b(setxor(1:length(b),nanid));
    s = s(setxor(1:length(s),nanid));
    line([0 0],[0 1],...
        'LineStyle','--',...
        'LineWidth',1.5,...
        'Color','k'), hold on
    errorbar(xs,b,s,'d','Color','k',...
        'LineWidth',1.5,...
        'markerfacecolor',[0,206,209]./255,...
        'markeredgecolor','k',...
        'markersize',12)
    axis([-80 100 0 0.5])
    set(gca,'LineWidth',1.5,'box','on',...
        'FontSize',14,'FontName','Arial',...
        'TickDir','out',...
        'TickLength',[0.009 0.009],...
        'YTick',0:0.1:0.5)
    ylabel('H_s (m)','FontSize',18)
    xlabel('Cross-shore distance (m)','FontSize',18)
%     export_fig([savefigdir 'XHs'],'-pdf','-nocrop')

end

