%Plot a timeseries of a|u|^3 and epsilon for a single deployment, single
%quadrat.

clear, close all
veldir = 'd:\Projects\Mekong_W2015\DataAnalysis\Paper2\AveragedVelocities\';
tkedir = 'd:\Projects\Mekong_W2015\DataAnalysis\Paper2\TKE\';
velfiles = dir(veldir);velfiles = {velfiles.name};
tkefiles = dir(tkedir);tkefiles = {tkefiles.name};
velfiles = velfiles(3:end);
tkefiles = tkefiles(3:end);

toprocess = [1 2 3 4 5 6 10];
for b = 1:length(toprocess)
    load([tkedir tkefiles{toprocess(b)}])
    load([veldir velfiles{toprocess(b)}])
    name = regexprep(velfiles{toprocess(b)},'Vave.mat','');
    if strfind(name,'HTA')
        name = ['FSS3_' name(end)];
    end
    %determine how many instruments are in the file
    fn = fieldnames(Avgs);
    [tp,~] = size(fn);
    if tp == 1
        g = 1;
    else
        g = 1:3;
    end
    for k = g                                                               %loop through instruments
        %%%Average TKE to the same length as Avgs
        n = length(Avgs.(fn{k}).time);
        beam1 = zeros(30,n);beam2 = zeros(30,n);beam3 = zeros(30,n);beam4 = zeros(30,n);
        for kk = 1:n-1
            td = find(Stat.(fn{k}).time >= Avgs.(fn{k}).time(kk) & Stat.(fn{k}).time <= Avgs.(fn{k}).time(kk+1));
            beam1(1:30,kk) = nanmean(Stat.(fn{k}).beam1.E(:,td),2);
            beam2(1:30,kk) = nanmean(Stat.(fn{k}).beam2.E(:,td),2);
            beam3(1:30,kk) = nanmean(Stat.(fn{k}).beam3.E(:,td),2);
            beam4(1:30,kk) = nanmean(Stat.(fn{k}).beam4.E(:,td),2);
        end
        %average beams together to get a single estimate of Eps per timestep
        Stat.(fn{k}).eps = (beam1+beam2+beam3+beam4)./4;
    end
    
    %load quadrat data
    quadir = 'd:\Projects\Documents\Writing\DataReports\ThirdAttempt\';
    order = [2 4 3 1]; %load files in order
    %%%%
    files = dir([quadir '*local40cm*.mat']);files = {files.name};
    for i = 1:length(files)
        load([quadir files{order(i)}])
        stage = regexp(files{order(i)},'.+_(.*).mat','tokens');
        fname = ['fou' char(stage{:})];
        vd.(fname).n = vegdat.n;
        vd.(fname).meanD = vegdat.meanD;
        vd.(fname).a = vegdat.a;
        vd.(fname).phi = vegdat.phi;
        vd.(fname).name = vegdat.Name;
    end
%     load([quadir 'Vegdat_1m.mat'])
%     vd.full.n = vegdat.n;
%     vd.full.meanD = vegdat.meanD;
%     vd.full.a = vegdat.a;
%     vd.full.phi = vegdat.phi;
%     vd.full.name = vegdat.Name;
    clear vegdat
    %Plot 3 instruments on the same figure (for F2F experiments, there will
    %only be two instruments since a = 0 for mudflat vectrinos!)
    % c = [68 245 255;204 20 137;255 219 0]./255;
    symb = {'^';'d';'s'};
    cmap = {'Purples';'Blues';'Greens'};
    pp = zeros(3,1);xmin = zeros(3,1);xmax = zeros(3,1);a=zeros(3,4);
    
    f1 = figure(b);
    set(f1,'PaperOrientation','portrait',...
        'position',[400 100   800 600]);
    set(gcf,'color','w','PaperPositionMode','auto','Renderer','zbuffer')
    
    vfn = fieldnames(vd);
    for i = 1:length(vfn)
        vid = strfind(vd.(vfn{i}).name,name);
        id = find(not(cellfun('isempty',vid)));
        a(:,i) = vd.(vfn{i}).a(id);
    end
    a = max(a,[],2); %find maximum a
    stage = vfn{1}(1:3);
    
    for i = g;
        ucubed = Avgs.(fn{i}).ucubed(7,:);
        eps = Stat.(fn{i}).eps(7,:);
        n = length(eps);
        c = brewermap(n,cmap{i});
        for ii = 1:n
            plot(a(i).*ucubed(ii),eps(ii),symb{i},...
                'MarkerFaceColor',c(ii,:),...
                'MarkerEdgeColor','k',...
                'MarkerSize',10); hold on
        end
        pp(i) = plot(0,0,symb{i},...
            'MarkerSize',10,...
            'MarkerFaceColor',[0.7 0.7 0.7],...
            'MarkerEdgeColor','k'); %for legend
        xmin(i) = min(a(i).*ucubed);
        xmax(i) = max(a(i).*ucubed);
    end
    xmin = min(xmin);xmax = max(xmax);
    set(gca,...
        'XScale','log',...
        'YScale','log',...
        'XLim',[xmin xmax],...
        'YLim',[1E-6 1E-1],...
        'TickDir','out',...
        'FontSize',14,...
        'FontName','Arial',...
        'LineWidth',1.5)
    c = brewermap(n/2,'Greys'); %for colorbar
    %averages every 30sec, divide n by 2 to get avgs in minutes
    colormap(c)
    h = colorbar('eastoutside');
    set(h,...
        'LineWidth',1.5,...
        'FontSize',14)
    ylabel(h,'Minutes Elapsed','FontSize',16)
    
    leg = legend(pp,fn);
    set(leg,'position',[0.7 0.2 0.05 0.05])
    xlabel('a|u|^3','FontSize',18)
    ylabel('\epsilon (Wkg^-^1)','FontSize',18)
    title(['Experiment: ' name])
    savefigdir = 'd:\Projects\Mekong_W2015\Figures\ShortFormat\';
    fname = [name '_aUcEps_all_' stage];
    export_fig([savefigdir fname],'-png','-nocrop');
end