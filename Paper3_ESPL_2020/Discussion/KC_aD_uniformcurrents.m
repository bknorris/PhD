%Plot KC number, Scour depth on one figure, and relationship between KC and
%ad (Nepf 1999, ML 2004) on another figure.
clear, close all
%define working paths
sfdir = 'd:\GradSchool\DataAnalysis\Paper3\Figures\';
dpath = '\DataAnalysis\Paper3\';
ypath = {'Mekong_F2014';'Mekong_W2015'};
%load run file to tell program which files & VPs to load
rdir = 'd:\GradSchool\DataAnalysis\Paper3\Analysis\';
fid = fopen([rdir 'ExpData_unicurrents.csv']);
rfile = textscan(fid,'%s%s%s%s%n%n%n%n%n%n%n%n%n%s','delimiter',',');
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
rTcr = rfile{12}; %Critical shear stress
rd50 = rfile{13}; %d50 of sediment
rinst = rfile{14}; %which adcp file to use
%initialize structure
KC = zeros(length(rexp),3);
SDobs = zeros(length(rexp),3);
SDuni1 = zeros(length(rexp),3);
SDuni2 = zeros(length(rexp),3);
SDuni3 = zeros(length(rexp),3);
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
            %load wave stats file
            wfile = dir([npath 'WaveStats\' folders{ii} '\', '*' rinst{vpid(j)} '_v2.mat']);
            load([npath 'WaveStats\' folders{ii} '\', wfile.name])
            %Calculate KC numbers
            ubr = data.(dfn{j}).ubr;
            tr = data.(dfn{j}).Tr;
            h = data.(dfn{j}).depth;
            tcr = rTcr(vpid(j));
            d50 = rd50(vpid(j));
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
            %KC
            KC(vpid(j),1) = nanmean(KCD);
            KC(vpid(j),2) = nanstd(KCD);
            KC(vpid(j),3) = nanstd(KCD);
            %SDobs
            Sd = abs(bdmin)./D;
            SDobs(vpid(j),1) = nanmean(Sd);
            SDobs(vpid(j),2) = nanstd(Sd);
            SDobs(vpid(j),3) = nanstd(Sd);
            %SDuni
            ys = nanmean(h)+(nanmean(abs(bdmin))); %max scour depth plus water depth
            rho = nanmean(wave.rho); %sea density
            uscrit = sqrt(tcr/rho); %critical shear velocity
            Uc = 5.75*uscrit*log10(5.53*(ys/d50)); %critical velocity, e.g., Melville & Coleman, 2000
            Fr = nanmean(ubr./sqrt(9.81*h)); %Froude number
            Frc = Uc/(sqrt(9.81*nanmean(h))); %critical Froude number
            %Try a bunch of different equations for estimating scour
            SDuni1(vpid(j),1) = 1.86*((nanmean(h)/d)^0.3)*abs(Fr-Frc)^0.25; %e.g., Jain & Fischer, 1980
            SDuni1(vpid(j),2) = 1.86*((nanmean(h)/d)^0.3)*abs(Fr-Frc)^0.25;
            SDuni1(vpid(j),3) = 1.86*((nanmean(h)/d)^0.3)*abs(Fr-Frc)^0.25;
            SDuni2(vpid(j),1) = 2.42*abs((2*nanmean(ubr)./Uc)-1)*((Uc^2)/(9.81*d))^(1/3); %e.g., Hancu, 1971
            SDuni2(vpid(j),2) = 2.42*abs((2*nanmean(ubr)./Uc)-1)*((Uc^2)/(9.81*d))^(1/3);
            SDuni2(vpid(j),3) = 2.42*abs((2*nanmean(ubr)./Uc)-1)*((Uc^2)/(9.81*d))^(1/3);
            SDuni3(vpid(j),1) = 1.84*((nanmean(h)/d)^0.3)*Frc^0.25; %e.g., Jain, 1981
            SDuni3(vpid(j),2) = 1.84*((nanmean(h)/d)^0.3)*Frc^0.25;
            SDuni3(vpid(j),3) = 1.84*((nanmean(h)/d)^0.3)*Frc^0.25;
        end
    end
end
KC = KC(~any(isnan(KC) | isinf(KC),2),:);
SDobs = SDobs(~any(isnan(SDobs) | isinf(SDobs),2),:);
SDuni1 = SDuni1(~any(isnan(SDuni1) | isinf(SDuni1),2),:);
SDuni2 = SDuni2(~any(isnan(SDuni2) | isinf(SDuni2),2),:);
SDuni3 = SDuni3(~any(isnan(SDuni3) | isinf(SDuni3),2),:);
rtide(strcmp(rarea,'mud')) = []; %remove mudflat values
rad(strcmp(rarea,'mud')) = [];
rarea(strcmp(rarea,'mud')) = [];
%Load the other data

% Construct figure
f1 = figure(1);
set(f1,'PaperOrientation','portrait',...
    'position',[400 100   800   400],...
    'renderer','painters');
p = zeros(3,3);
p(1,:) = plot(SDobs,SDuni1,'ok','markerfacecolor','r');hold on
p(2,:) = plot(SDobs,SDuni2,'^k','markerfacecolor','b');
p(3,:) = plot(SDobs,SDuni3,'sk','markerfacecolor','g');
leg = legend(p(:,1),{'Jain & Fischer, 1980';'Hancu, 1971';'Jain, 1981'});
% set(gca,...
%     'xlim',[4 150],...
%     'ylim',[0.005 2])
%Global Adjustments
set(leg,'position',[0.53 0.205 0.05 0.05])
xlabel('S/D_{obs} [-]')
ylabel('S/D_{uni} [-]')
prettyfigures('text',12,'labels',13,'box',1,'gcolor','k')
export_fig([sfdir 'SDobs_SDuni'],'-pdf','-nocrop')

