%estimate ripple properties based on hydrodynamics, Coco et al. 2007

clear, close all
%Define working paths
dpath = '\DataAnalysis\Paper3\';
ypath = {'Mekong_F2014';'Mekong_W2015'};
%load run file to tell program which files & VPs to load
rdir = 'g:\GradSchool\DataAnalysis\Paper3\ExperimentalDesign\';
fid = fopen([rdir 'RippleRunfile.csv']);
rfile = textscan(fid,'%s%s%s%n%n%n%n%n%s','delimiter',',');
rexp = rfile{1}; %experiment name
rdate = rfile{2}; %date of experiment
rfld = rfile{3}; %area of forest
rd50 = rfile{4}; %sediment d50 (m)
rd90 = rfile{5}; %sediment d90 (m)
rsand = rfile{6}; %percent sand in sample (%)
rmud = rfile{7}; %percent mud in sample (%)
rcrit = rfile{8}; %critical shear stress from Critical_sed_params.m
aqfil = rfile{9}; %which aquadopp to load
%initialize data file
dat.mud = [];
dat.fringe = [];
dat.forest = [];
%% Loop through folder structure
for i = 2
	npath = ['f:\' ypath{i} dpath];
    folders = dir([npath '\VPs\']);
    folders = {folders(3:end).name};
    for ii = 1:length(folders)
        if any(strcmp(folders{ii},rdate))
            vpid = strcmp(folders{ii},rdate);vpid = find(vpid);
            disp(['Loading files from ' npath 'WaveStats\' folders{ii} '\'])
            afile = dir([npath 'WaveStats\' folders{ii} '\','*.mat']);afile = {afile.name};
        else
            continue
        end
        for j = 1:length(vpid)
            %% Load files
            disp(['Processing ' rfld{vpid(j)}])
            aqid = strfind(afile, aqfil{vpid(j)});aqid = find(~cellfun(@isempty,aqid));
            load([npath 'WaveStats\' folders{ii} '\' afile{aqid}])
            
            %% Define parameters
            
            Aw = wave.Uwrms./wave.omegar; %rms wave orbital excursion
            Uw = wave.Uwrms; %rms wave orbital velocity
            nu = 1.05E-6; %viscosity of seawater
            g = 9.81; 
            
            %sediment properties
            d50 = rd50(vpid(j));
            pmapos = 0.13+(0.21/(d50*1000+0.002)^0.21); %eq 4 from Wu & Wang, 2006
            rhos = 2.650; %from Coco et al. 2007
%             rhod = (1-pmapos)*rhom; %dry mixture density based off of d50.
            rhow = 1.025; %density of seawater
            
            %Coco et al. 2007 method
            eta = zeros(length(Uw),1);
            lambda = zeros(length(Uw),1);
            for k = 1:length(Uw)
                %eq 23 from Coco et al. 2007
                X = (4*nu*Uw(k)^2)/(d50*(((rhos/rhow)-1)*g*d50)^(1.5));
                if X <= 2
                    eta(k) = (0.30*(X^-0.39));
                    lambda(k) = (1.96*(X^-0.28));
                elseif X >= 2
                    eta(k) = (0.45*(X^-0.99));
                    lambda(k) = (2.71*(X^-0.75));
                end
            end
            
            %calculate bedload flux per unit width and bedform migration
            %velocity (e.g., Coco et al. 2007, Smyth & Li, 2005)
            fw = (2*wave.tauw)./(wave.rho.*(wave.ubr.^2)); %get fw from wave data
            theta = (0.5*fw.*(Uw.^2))./((rhos-rhow)*g*d50*(1-((pi.*eta')./lambda')).^2); %shields parameter
            thetac = rcrit(vpid(j))/(rhow*((rhos/rhow)-1)*g*d50*1000);
            %Dstar = d50*(((g*((rhos-rhow)/rhow))/nu^2)^(1/3)); %
            qb = 9.1*((g*((rhos-rhow)/rhow)*(1000*d50^3))^(1/2))*(abs(theta-thetac).^1.78); %sed flux per unit width
            Um = qb./(eta'*(1-pmapos)); %bedform migration velocity
            
            %% Save to structure
            if isempty(dat.(rfld{vpid(j)}))
                dat.(rfld{vpid(j)}).time = [];
                dat.(rfld{vpid(j)}).etA = [];
                dat.(rfld{vpid(j)}).lambdA = [];
                dat.(rfld{vpid(j)}).eta = [];
                dat.(rfld{vpid(j)}).lambda = [];
                dat.(rfld{vpid(j)}).Hs = [];
                dat.(rfld{vpid(j)}).Um = [];
                dat.(rfld{vpid(j)}).Uw = [];
            else
                dat.(rfld{vpid(j)}).time = [dat.(rfld{vpid(j)}).time; wave.time2']; %time
                dat.(rfld{vpid(j)}).etA = [dat.(rfld{vpid(j)}).etA; eta]; %eta/Aw
                dat.(rfld{vpid(j)}).lambdA = [dat.(rfld{vpid(j)}).lambdA; lambda]; %lambda/Aw
                dat.(rfld{vpid(j)}).eta = [dat.(rfld{vpid(j)}).eta; eta.*Aw']; %eta
                dat.(rfld{vpid(j)}).lambda = [dat.(rfld{vpid(j)}).lambda; lambda.*Aw']; %lambda
                dat.(rfld{vpid(j)}).Hs = [dat.(rfld{vpid(j)}).Hs; wave.Hs']; %sig. wave height
                dat.(rfld{vpid(j)}).Um = [dat.(rfld{vpid(j)}).Um; Um']; %bedform velocity
                dat.(rfld{vpid(j)}).Uw = [dat.(rfld{vpid(j)}).Uw; Uw']; %wave orbital velocity
            end
        end
    end
end

figure(1)
sp(1) = subplot(311);
hist(dat.mud.lambda,10),title('Mudflat')
sp(2) = subplot(312);
hist(dat.fringe.lambda,10),title('Fringe')
sp(3) = subplot(313);
hist(dat.forest.lambda,10),title('Forest')
xlabel('\lambda [m]')
set(sp,'xlim',[0 0.14])
% 
figure(2)
sp(1) = subplot(311);
hist(dat.mud.eta,10),title('Mudflat')
sp(2) = subplot(312);
hist(dat.fringe.eta,10),title('Fringe')
sp(3) = subplot(313);
hist(dat.forest.eta,10),title('Forest')
xlabel('\eta [m]')
% set(sp,'xlim',[0 4E-4],'ylim',[0 6])
figure(3)
sp(1) = subplot(311);
hist(dat.mud.Uw,10),title('Mudflat')
sp(2) = subplot(312);
hist(dat.fringe.Uw,10),title('Fringe')
sp(3) = subplot(313);
hist(dat.forest.Uw,10),title('Forest')
xlabel('U_w [m/s]')
set(sp,'xlim',[0 0.3])

figure(4)
sp(1) = subplot(311);
title('Significant Wave Height')
plot(dat.mud.time,dat.mud.Hs,'o')
hold on
plot(dat.fringe.time,dat.fringe.Hs,'^r')
plot(dat.forest.time,dat.forest.Hs,'sg')

sp(2) = subplot(312);
title('Ripple Wavelength')
plot(dat.mud.time,dat.mud.lambda,'o')
hold on
plot(dat.fringe.time,dat.fringe.lambda,'^r')
plot(dat.forest.time,dat.forest.lambda,'sg')

sp(3) = subplot(313);
title('Ripple Migration Velocity')
plot(dat.mud.time,dat.mud.Um*6000,'o')
hold on
plot(dat.fringe.time,dat.fringe.Um*6000,'^r')
plot(dat.forest.time,dat.forest.Um*6000,'sg')
set([sp(1) sp(2)],'xticklabel',[])
ylabel('cm/hr')
datetickzoom('x','dd-mm HH:MM','keepticks','keeplimits')
linkaxes(sp,'x')
