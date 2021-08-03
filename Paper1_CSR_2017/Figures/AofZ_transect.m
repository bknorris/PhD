%Plot a (frontal area density) as a function of z for every x-shore
%position along the transect line; Paper 1. This task is actually fairly
%complicated; I need to load in each of the quad reconstructions and
%calculate a for each layer.
clear
%If the file exits all ready, skip the loading step
datdir = 'd:\Projects\Mekong_W2015\DataAnalysis\DataReports\Vegetation\';
if exist([datdir 'VegDat_by_height_1m.mat'],'file') > 0
    load([datdir 'VegDat_by_height_1m.mat'])
else
    dirs{1} = 'd:\Projects\Mekong_W2015\DataAnalysis\DataReports\Vegetation\2014\';
    dirs{2} = 'd:\Projects\Mekong_W2015\DataAnalysis\DataReports\Vegetation\2015\';
    yr = {'four';'five'};
    files = cell(2,10);
    dir1 = dir(dirs{1});files(1,1:9) = {dir1.name};
    dir2 = dir(dirs{2});files(2,1:10) = {dir2.name};
    files(:,1:2) = []; %remove up/down directories
    for i = 1:2
        for ii = 1:length(files)
            if ~isempty(files{i,ii})
                disp(['Loading: ' dirs{i} files{i,ii}])
                load([dirs{i} files{i,ii}])
                name = regexprep(files{i,ii},'stats_fixed.mat','');
                DATA = rmfield(DATA,'info');
                fn = fieldnames(DATA);
                v = length(fn);
                h = 0.005;
                z = h:h:v*h;
                
                [vv,~] = size(DATA.layer1.Rcenter);
                Acir = zeros(v,vv);acir = zeros(v,1);
                D = zeros(v,vv);d = zeros(v,1);
                n = zeros(v,1);
                for j = 1:v
                    %calculate veg geometries
                    [ll,~] = size(DATA.(fn{j}).Rradii);
                    for jj = 1:ll
                        r = mean(DATA.(fn{j}).Rradii(jj,:));
                        D(j,jj) = r*2;
                        Acir(j,jj) = pi*r^2;
                    end
                    n(j) = length(DATA.(fn{j}).Rradii);
                    d(j) = mean(D(j,:));
                    acir(j) = sum(Acir(j,:));
                end
                area = 1; %m^2
                a = (max(n)/area).*d;
                phi = acir/area;
                
                %save to structure
                veg.(yr{i}).(name).z = z;
                veg.(yr{i}).(name).n = n;
                veg.(yr{i}).(name).d = d;
                veg.(yr{i}).(name).a = a;
                veg.(yr{i}).(name).phi = phi;
            else
                continue
            end
        end
    end
    save([datdir 'VegDat_by_height_1m'],'veg','-v7.3')
end
%Get VP positions along x-shore transect, then plot the figure
load('D:\Projects\Mekong_W2015\DataAnalysis\DataReports\Paper1_QuadAnalysis\ThirdAttempt\Vegdat_1m.mat')
f1 = figure(1);
set(f1,'PaperOrientation','portrait',...
    'position',[400 100   1000 300]);
set(gcf,'color','w','PaperPositionMode','auto')
fn = fieldnames(veg);
c = 1;
cmap = [238,114,137;241,116,94;226,147,110;233,126,57;224,158,55;...
    187,158,79;234,211,125;188,177,54;215,230,102;133,190,69;...
    182,225,147;103,187,101;109,170,109;102,226,128;89,225,168;...
    83,219,197;62,198,227;95,158,234;159,143,238]./255;
p = zeros(19,1);
%plot fringe delineations
xfringe = [-10;20]*ones(1,10);
yfringe = linspace(0,1,10);
plot(xfringe,yfringe,'--k','linewidth',1.5),hold on
%plot Mudflat quads
Xpos = [-37.17 -60.96 -50.54 -50.35];
for i = 1:length(Xpos)
    scale = 0.1; %scale adjustment factor for plot (increment between ticks)
    a = zeros(10,1);
    z = linspace(0,1,10);
    x = Xpos(i)+(a/scale);
%     p(c) = plot(x,z,'color',cmap(c,:),'linewidth',1.5);
    p(c) = plot(Xpos(i),0.5,'o','color',cmap(c,:),'markerfacecolor',cmap(c,:));
    text(Xpos(i),0.51,num2str(c))
    c = c+1;
end
%Plot Fringe-Forest quads
for i = 1:2
    dn = fieldnames(veg.(fn{i}));
    for ii = 1:length(dn)
        if i == 1
            txt = 'first'; %there are 2 Q4's
        else
            txt = 'last';
        end
        Qid = find(strcmp(dn{ii},vegdat.Qname),1,txt);
        Xpos = vegdat.Xshore(Qid);
        a = veg.(fn{i}).(dn{ii}).a;
        z = veg.(fn{i}).(dn{ii}).z;
        x = Xpos+(a/scale);
        p(c) = plot(x,z,'color',cmap(c,:),'linewidth',1.5);
%         p(c) = plot(Xpos,0.5,'o','color',cmap(c,:),'markerfacecolor',cmap(c,:));
%         text(Xpos,0.51,num2str(c))
        c = c+1;
    end
end
%Plot a scale bar
plot([82 82+(1/scale)],[0.65 0.65],'k','Linewidth',2)
% qname = {'Mudflat';'Mudflat';'Mudflat';'Mudflat';...
%     'Q1';'Q2';'Q3a';'Q3b';'Q3c';'Q4';'Q5';'Q6';...
%     'Q7';'Q8';'Q10';'Q11';'Q12';'Q13'};
% leg = legend(p,qname);
set(gca,...
    'xlim',[-80 100],...
    'ylim',[0 0.8],...
    'linewidth',1.5,...
    'box','on',...
    'FontSize',14,...
    'FontName','Arial',...
    'TickDir','out')
xlabel('Cross Shore Distance (m)','fontsize',18)
ylabel('z (m)','fontsize',18)
figdir = 'd:\Projects\Mekong_W2015\Figures\Paper1\';
% export_fig([figdir 'AofZ_transect'],'-pdf','-nocrop')

