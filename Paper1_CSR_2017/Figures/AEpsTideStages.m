clear
%load the vegdat files and name based on the size and stage of tide
vpdir = 'd:\Projects\Documents\Writing\DataReports\SecondAttempt\';
vpdir2 = 'd:\Projects\Documents\Writing\DataReports\ThirdAttempt\';
vd = struct(); %big vegdat structure
order = [2 4 3 1]; %load files in order

%%%10cm%%%
files = dir([vpdir '*local10cm*.mat']);files = {files.name};
for i = 1:length(files)
    load([vpdir files{order(i)}])
    stage = regexp(files{order(i)},'.+_(.*)_vol.mat','tokens');
    fname = ['ten' char(stage{:})];
    vd.(fname).n = vegdat.n;
    vd.(fname).meanD = vegdat.meanD;
    vd.(fname).a = vegdat.a_vol;
    vd.(fname).phi = vegdat.phi_vol;
end
%%%20cm%%%
files = dir([vpdir2 '*local20cm*.mat']);files = {files.name};
for i = 1:length(files)
    load([vpdir2 files{order(i)}])
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
    stage = regexp(files{order(i)},'.+_(.*)_vol.mat','tokens');
    fname = ['thr' char(stage{:})];
    vd.(fname).n = vegdat.n;
    vd.(fname).meanD = vegdat.meanD;
    vd.(fname).a = vegdat.a_vol;
    vd.(fname).phi = vegdat.phi_vol;
end
%%%40cm%%%
files = dir([vpdir '*local40cm*.mat']);files = {files.name};
for i = 1:length(files)
    load([vpdir files{order(i)}])
    stage = regexp(files{order(i)},'.+_(.*)_vol.mat','tokens');
    fname = ['fou' char(stage{:})];
    vd.(fname).n = vegdat.n;
    vd.(fname).meanD = vegdat.meanD;
    vd.(fname).a = vegdat.a_vol;
    vd.(fname).phi = vegdat.phi_vol;
end
%%%50cm%%%
files = dir([vpdir '*local50cm*.mat']);files = {files.name};
for i = 1:length(files)
    load([vpdir files{order(i)}])
    stage = regexp(files{order(i)},'.+_(.*)_vol.mat','tokens');
    fname = ['fif' char(stage{:})];
    vd.(fname).n = vegdat.n;
    vd.(fname).meanD = vegdat.meanD;
    vd.(fname).a = vegdat.a_vol;
    vd.(fname).phi = vegdat.phi_vol;
end
%%%1m%%%
load([vpdir 'Vegdat_local1m_vol.mat'])
vd.full.n = vegdat.n;
vd.full.meanD = vegdat.meanD;
vd.full.a = vegdat.a_vol;
vd.full.phi = vegdat.phi_vol;
vd.X = vegdat.Xshore;
clear vegdat
vegdat = vd;

%%%Velocity Statistics%%%
CatTurbData_v2
fn = fieldnames(vegdat);
E = [veldat.four.E; veldat.five.E];
std = [veldat.four.Estd; veldat.five.Estd];
ucubed = [veldat.four.ucubed; veldat.five.ucubed];
H = [veldat.four.depth; veldat.five.depth];
Hs = [veldat.four.wrms; veldat.five.wrms].*sqrt(2);
gamma = Hs./H;
ucubed(ucubed > 0.005) = 0.002212;
savefigdir = 'd:\Projects\Mekong_W2015\Figures\Paper2\Phi-Eps\';

plota = 0;
plotaucubed = 1;
break
if plota == 1
    %%%X and A for a single scale, including LL to HH (20cm^2)
    symb = {'^','o','d','s'};
    lins = {':','-.','--','-'};
    f2 = figure(2);
    set(f2,'PaperOrientation','portrait',...
        'position',[400 100   1000 800]);
    set(gcf,'color','w','PaperPositionMode','auto')
    sp(1) = subplot(311);
    vid = strfind(fn,'twe');
    id = find(not(cellfun('isempty',vid)));
    step = 10;
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
    A = zeros(54,4);
    for i = 1:length(id)
        A(:,i) = vegdat.(fn{id(i)}).a;
    end
    A = nanmean(A,2);
    X = vegdat.X;
    xs = min(X):step:max(X);
    [b,~,s] = bindata(X,A,xs);
    nanid = find(isnan(b));
    %remove NaNs
    xs = xs(setxor(1:length(xs),nanid));
    b = b(setxor(1:length(b),nanid));
    s = s(setxor(1:length(s),nanid));
    s(5)= 0.016;
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
        'YLim',[-0.002 0.05],...
        'YTick',0:0.01:0.05,...
        'TickDir','out')
    ylabel('a','FontSize',18)
    
    sp(3) = subplot(313);
    c = brewermap(4,'YlGnBu');colormap(c)
    c(1,:) = [255 235 0]./255;
    p = zeros(1,4);
    for i = 1:4
        [b,~,s] = bindata(X,E(:,i),xs);
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
    end
    set(gca,'LineWidth',1.5,...
        'FontSize',14,...
        'FontName','Arial',...
        'YLim',[-1E-4 4E-3],...
        'TickDir','out')
    leg = legend(p,{'Low';'Mid-Low';'Mid-High';'High'});
    set(leg,'position',[0.88 0.25 0.05 0.05],...
        'LineWidth',1.5)
    ylabel('\epsilon (Wkg^-^1)','FontSize',18)
    xlabel('Cross-shore distance (m)','FontSize',18)
    
    set(sp(1),'Position',[0.1 0.68 0.84 0.23],...
        'XTickLabel',[])
    set(sp(2),'Position',[0.1 0.4 0.84 0.23],...
        'XTickLabel',[])
    set(sp(3),'Position',[0.1 0.12 0.84 0.23])
    hold off
    %     export_fig([savefigdir 'Xa_Xeps_LL_HH'],'-pdf','-nocrop')
    
    %%%Subdivide the forest into two sections, plot a/Eps for both
    %%%sections%%%
    f3 = figure(3);
    set(f3,'PaperOrientation','portrait',...
        'position',[400 100   800 600]);
    set(gcf,'color','w','PaperPositionMode','auto')
    c = brewermap(4,'YlGnBu');colormap(c)
    c(1,:) = [255 235 0]./255;
    vid = strfind(fn,'twe');
    id = find(not(cellfun('isempty',vid)));
    step = [0.007 0.008 0.009 0.01];
    p = zeros(1,4);
    for i = 1:length(id)
        A = vegdat.(fn{id(i)}).a;
        zid = find(A == 0);
        zid = setxor(1:54,zid);
        Ez = E(zid);Ez(Ez > 1E-3) = NaN;
        az = A(zid);az(az > 0.05) = NaN;
        xs = min(az):step(i):max(az);
        [b,~,s] = bindata(az,Ez,xs);
        nanid = find(isnan(b));
        %remove NaNs
        xs = xs(setxor(1:length(xs),nanid));
        b = b(setxor(1:length(b),nanid));
        s = s(setxor(1:length(s),nanid));                %edit 15/07/2016: plot all points (for Julia)
        plot(xs,b,lins{i},'LineWidth',1.5,...
            'Color',c(i,:)), hold on
        %         q(i) = plot(phiz,Ez,symb{i},...
        %             'Markersize',6,...
        %             'MarkerFaceColor',[0.5 0.5 0.5],...
        %             'MarkerEdgeColor','k');
        p(i) = errorbar(xs,b,s,symb{i},...
            'Color','k',...
            'LineWidth',1.5,...
            'markerfacecolor',c(i,:),...
            'markeredgecolor','k',...
            'markersize',10);
        
        hold on
    end
    set(gca,'LineWidth',1.5,...
        'FontSize',13,...
        'FontName','Arial',...
        'YLim',[-1E-4 1.1E-3],...
        'XLim',[-1E-3 0.05],...
        'XTick',0:0.01:0.05,...
        'YTick',0:1E-4:1E-3,...
        'TickDir','out')
    leg = legend(p,{'Low';'Mid-Low';'Mid-High';'High'});
    set(leg,'position',[0.22 0.78 0.05 0.05],...
        'LineWidth',1.5,'box','off')
    ylabel('\epsilon (Wkg^-^1)','FontSize',18)
    xlabel('a','FontSize',18)
    %     export_fig([savefigdir '20cm_example_allsamples_a'],'-pdf','-nocrop')
    
        f1 = figure(1);
    set(f1,'PaperOrientation','portrait',...
        'position',[400 100   600 500]);
        set(gcf,'color','w','PaperPositionMode','auto')

    c = brewermap(4,'YlGnBu');colormap(c)
    c(1,:) = [255 235 0]./255;
    vid = strfind(fn,'twe');
    t = {'Low';'Mid-Low';'Mid-High';'High'};
    txt = [0.018 0.014 0.01 0.004];
    id = find(not(cellfun('isempty',vid)));
    p = zeros(1,4);
    for i = 1:length(id)
        A = vegdat.(fn{id(i)}).a;
        A(A == 0) = NaN;
        nid = ~isnan(A);
        az = A(nid);
        phi = vegdat.(fn{id(i)}).phi;
        phiz = phi(nid);
        %linreg
        pf = polyfit(az,phiz,1);
        stats = regstats(phiz,az,'linear','fstat');
        disp([t{i} ' line slope ' num2str(pf(1))])
        disp([t{i} ' p-value ' num2str(stats.fstat.pval)])
        pv = polyval(pf,az);
        resid = phiz-pv;
        SSresid = sum(resid.^2);
        SStotal = (length(phiz)-1)*var(phiz);
        rsq = 1-SSresid/SStotal;
        p(i) = plot(az,phiz,'o',...
            'LineWidth',1.5,...
            'markerfacecolor',c(i,:),...
            'markeredgecolor','k',...
            'markersize',10);
        hold on
        plot(az,pv,'LineWidth',1.5,...
            'Color',c(i,:))
        text(4,txt(i),['R-sqd, ' t{i} ': ' sprintf('%.2f',rsq)])
    end
    set(gca,'LineWidth',1.5,...
        'FontSize',14,...
        'FontName','Arial',...
        'TickDir','out')
    xlabel('a (m^-^1)','FontSize',18)
    ylabel('\epsilon (Wkg^-^1)','FontSize',18)
%     export_fig([savefigdir 'a_phi'],'-jpeg','-nocrop')

end
if plotaucubed == 1
    %%%Subdivide the forest into two sections, plot a/Eps for both
    %%%sections%%%
    symb = {'^','o','d','s'};
    lins = {':','-.','--','-'};
    f1 = figure(1);
    set(f1,'PaperOrientation','portrait',...
        'position',[400 100   800 600]);
    set(gcf,'color','w','PaperPositionMode','auto')
    c = brewermap(4,'YlGnBu');colormap(c)
    c(1,:) = [255 235 0]./255;
    vid = strfind(fn,'twe');
    id = find(not(cellfun('isempty',vid)));
    step = [1E-6 1E-6 3E-7 2E-7];
    p = zeros(1,4);
    for i = 1:length(id)
        A = vegdat.(fn{id(i)}).a.*ucubed(:,i);
        zid = find(A == 0);
        zid = setxor(1:54,zid);
        Ez = E(zid);Ez(Ez > 1E-3) = NaN;
        az = A(zid);
        xs = min(az):step(i):max(az);
        [b,~,s] = bindata(az,Ez,xs);
        nanid = find(isnan(b));
        %remove NaNs
        xs = xs(setxor(1:length(xs),nanid));
        b = b(setxor(1:length(b),nanid));
        s = s(setxor(1:length(s),nanid));                %edit 15/07/2016: plot all points (for Julia)
        plot(xs,b,lins{i},'LineWidth',1.5,...
            'Color',c(i,:)), hold on
        %         q(i) = plot(phiz,Ez,symb{i},...
        %             'Markersize',6,...
        %             'MarkerFaceColor',[0.5 0.5 0.5],...
        %             'MarkerEdgeColor','k');
        p(i) = errorbar(xs,b,s,symb{i},...
            'Color','k',...
            'LineWidth',1.5,...
            'markerfacecolor',c(i,:),...
            'markeredgecolor','k',...
            'markersize',10);
        %         plot(az,Ez,'o','color',c(i,:))
        
        hold on
    end
    set(gca,'LineWidth',1.5,...
        'FontSize',13,...
        'FontName','Arial',...
        'YLim',[-2E-5 6E-4],...
        'XLim',[-5E-7 5E-6],...
        'TickDir','out')
    leg = legend(p,{'Low';'Mid-Low';'Mid-High';'High'});
    set(leg,'position',[0.72 0.78 0.05 0.05],...
        'LineWidth',1.5,'box','off')
    ylabel('\epsilon (Wkg^-^1)','FontSize',18)
    xlabel('a*|u^3| (m^2s^-^3)','FontSize',18)
    export_fig([savefigdir '20cm_auc_zoom'],'-pdf','-nocrop')
    
    %plot a*ucubed by # of samples
    f2 = figure(2);
    set(f2,'PaperOrientation','portrait',...
        'position',[400 100   600 500]);
    set(gcf,'color','w','PaperPositionMode','auto')
    c = brewermap(4,'YlGnBu');colormap(c)
    c(1,:) = [255 235 0]./255;
    vid = strfind(fn,'twe');
    id = find(not(cellfun('isempty',vid)));
    p = zeros(1,4);
    for i = 1:length(id)
        A = vegdat.(fn{id(i)}).a.*ucubed(:,i);
        p(i) = plot(1:54,A,'-o',...
            'LineWidth',1.5,...
            'Color',c(i,:),...
            'MarkerEdgeColor','k',...
            'MarkerFaceColor',c(i,:));
        hold on
    end
    set(gca,'LineWidth',1.5,...
        'FontSize',14,...
        'FontName','Arial',...
        'YLim',[-1E-6 8E-5],...
        'TickDir','out')
    leg = legend(p,{'Low';'Mid-Low';'Mid-High';'High'});
       set(leg,'position',[0.22 0.78 0.05 0.05],...
        'LineWidth',1.5,'box','off')
        xlabel('Quadrat #','FontSize',18)
    ylabel('a*|u^3| (m^2s^-^3)','FontSize',18)
    export_fig([savefigdir 'Auc_byquad'],'-jpeg','-nocrop')
end
