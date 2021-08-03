%Plot KC number, Scour depth on one figure, and relationship between KC and
%ad (Nepf 1999, ML 2004) on another figure.
clear, close all
%define working paths
sfdir = 'd:\MemoryStick\GradSchool\DataAnalysis\Paper3\Figures\';
dpath = '\DataAnalysis\Paper3\';
ypath = {'Mekong_F2014';'Mekong_W2015'};
%load run file to tell program which files & VPs to load
rdir = 'd:\MemoryStick\GradSchool\DataAnalysis\Paper3\Analysis\';
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
    npath = ['d:\' ypath{i} dpath];
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
            KCD = (ubr.*tr);
            %KCD
            KC(vpid(j),1) = nanmean(KCD);
            KC(vpid(j),2) = nanstd(KCD);
            KC(vpid(j),3) = nanstd(KCD);
            %Sd
            Sd = abs(bdmin);
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
%Load the other data
fid = fopen([rdir 'ML2004_Fig4.csv']);
rfile = textscan(fid,'%n%n','delimiter',',');
mlKC = rfile{1};mlCd = rfile{2};
fid = fopen([rdir 'Nepf1999_Fig6.csv']);
rfile = textscan(fid,'%n%n','delimiter',',');
nfad = rfile{1};nfCd = rfile{2};
%Cd = 0.47exp(-0.052KC) -> from M&L_2004
%Cd = -0.2213log(ad)+0.0434 -> from Nepf 1999
KCexp = log(0.47./nfCd)./0.052;
ADexp = log(-0.365./(0.261-mlCd))./12.7;
xs = [real(ADexp); rad; nfad(KCexp>0)];
ys = [mlKC; KC(:,1); KCexp(KCexp>0)];

% Construct figure
f1 = figure(1);
set(f1,'PaperOrientation','portrait',...
    'position',[400 100   400   400],...
    'renderer','painters');
symb = {'^';'s'};
line = {'-.';'--'};
cl = [0.9 0.9 0.9;0.4 0.4 0.4];
pi = zeros(2,1);
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
    'xlim',[0.1 2],...
    'ylim',[0.00001 0.1])
%Global Adjustments
set(leg,'position',[0.75 0.205 0.05 0.05])
xlabel('U*T [-]')
ylabel('S [-]')
prettyfigures('text',12,'labels',13,'box',1,'gcolor','k')
set(f1,'units','inches');
pos = get(f1,'Position');
set(f1,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
% print(f1,[sfdir 'KC_SD_ad'],'-dpdf','-r0')