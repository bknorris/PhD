%Plot KC number, Scour depth on one figure, and relationship between KC and
%ad (Nepf 1999, ML 2004) on another figure.
clear, close all
%define working paths
sfdir = 'e:\MemoryStick\GradSchool\DataAnalysis\Paper3\Figures\';
dpath = '\DataAnalysis\Paper3\';
ypath = {'Mekong_F2014';'Mekong_W2015'};
%load run file to tell program which files & VPs to load
rdir = 'e:\MemoryStick\GradSchool\DataAnalysis\Paper3\Analysis\';
fid = fopen([rdir 'ExperimentalData.csv']);
rfile = textscan(fid,'%s%s%s%s%n%n%n%n%n%n%n','delimiter',',');
rexp = rfile{1}; %exp name
rdate = rfile{2}; %date of exp
rarea = rfile{3}; %area of exp
rtide = rfile{4}; %tidal stage
rhab = rfile{5}; %vectrino hab
rn = rfile{6}; %num stems
rd = rfile{7}; %stem diameter (d)
rD = rfile{8}; %patch diameter (D)
ra = rfile{9}; %frontal area density
rad = rfile{10}; %ad
raD = rfile{11}; %aD
%initialize structure
KC = zeros(length(rexp),3);
SD = zeros(length(rexp),3);
for i = 1:2
    npath = ['e:\' ypath{i} dpath];
    folders = dir([npath '\CmbData\']);
    folders = {folders(3:end).name};
    for ii = 1:length(folders)
        if any(strcmp(folders{ii},rdate))
            vpid = strcmp(folders{ii},rdate);vpid = find(vpid);
            load([npath 'CmbData\' folders{ii} '\','AllEventData.mat']);
        else
            continue
        end
        dfn = fieldnames(data);
        for j = 1:length(dfn)
            %Calculate KC numbers
            ubr = data.(dfn{j}).ubr;
            tr = data.(dfn{j}).Tr;
            bd = data.(dfn{j}).bdist;
            bst = find(~isnan(data.(dfn{j}).bdist),1,'first');
            bd = bd-rhab(vpid(j));
            %find min 500 values
            bdmin = zeros(500,1);
            for jj = 1:500
                [bdmin(jj),idx] = min(bd);
                % remove for the next iteration the last smallest value:
                bd(idx) = [];
            end
            d = rd(vpid(j));
            D = rD(vpid(j));
            KCD = (ubr.*tr)/d;
            %KCD
            KC(vpid(j),1) = nanmean(KCD);
            KC(vpid(j),2) = nanstd(KCD);
            KC(vpid(j),3) = nanstd(KCD);
            %Sd
            Sd = abs(bdmin)./D;
            SD(vpid(j),1) = nanmean(Sd);
            SD(vpid(j),2) = nanstd(Sd);
            SD(vpid(j),3) = nanstd(Sd);
        end
    end
end
KC = KC(~any(isnan(KC) | isinf(KC),2),:);
SD = SD(~any(isnan(SD) | isinf(SD),2),:);
rtide(strcmp(rarea,'mud')) = []; %remove mudflat values
rad(strcmp(rarea,'mud')) = [];
rn(strcmp(rarea,'mud')) = [];

% Construct figure
f1 = figure(1);
set(f1,'PaperOrientation','portrait',...
    'position',[400 100   800   400],...
    'renderer','painters');
symb = {'^';'s'};
line = {'-.';'--'};
cl = [0.9 0.9 0.9;0.4 0.4 0.4];
pi = zeros(2,1);
sp(1) = subplot(121);
xx = linspace(6,150,1000);
sfkc = 1.3*(1-exp(-0.03*(xx-6)));
sf = loglog(xx,sfkc,'-',...
    'color','k',...
    'linewidth',1.5);hold on
stage = {'flood','ebb'};
for i = 1:2
    idx = strcmp(rtide,stage{i});
    eb = ploterr(KC(idx,1),SD(idx,1),[KC(idx,2) KC(idx,3)],[SD(idx,2) SD(idx,3)],...
        symb{i},'logxy','abshhxy',0.04);hold on
    set(eb,...
        'color','k',...
        'linewidth',1,...
        'markerfacecolor',cl(i,:),...
        'markeredgecolor','k',...
        'markersize',6);
    pi(i) = eb(1);
end
leg = legend([pi; sf],{'Flood';'Ebb';'S1992'});
set(gca,...
    'xlim',[4 150],...
    'ylim',[0.01 5])
sp(2) = subplot(122);
hold on
%Plot datasets
for i = 1:2
    idx = strcmp(rtide,stage{i});
    eb = ploterr(rn(idx,1),SD(idx,1),[zeros(nnz(idx),1) zeros(nnz(idx),1)],[SD(idx,2) SD(idx,3)],...
        symb{i},'logxy','abshhxy',0.04);
    set(eb,...
        'color','k',...
        'linewidth',1,...
        'markerfacecolor',cl(i,:),...
        'markeredgecolor','k',...
        'markersize',6);
    pi(i) = eb(1);
end
leg2 = legend(pi,{'Flood';'Ebb'});
set(gca,...
    'xlim',[0 180],...
    'ylim',[0.01 0.7])
%Global Adjustments
set(sp(1),'position',[0.12 0.15 0.35 0.8])
set(sp(2),'position',[0.58 0.15 0.35 0.8])
set(leg,'position',[0.38 0.205 0.05 0.05])
set(leg2,'position',[0.84 0.85 0.05 0.05])
xlabel(sp(1),'KC [-]')
ylabel(sp(1),'S/D [-]')
xlabel(sp(2),'n [m^{-2}]')
ylabel(sp(2),'S/D [-]')
prettyfigures('text',12,'labels',13,'box',1,'gcolor','k')
set(f1,'units','inches');
pos = get(f1,'Position');
set(f1,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(f1,[sfdir 'KC_SD_n_SD'],'-dpdf','-r0')