%plot bottom trace, turbulent dissipation, bed shear stress, velocity
%squared
clear, close all
dpath = 'd:\Projects\Mekong_W2015\DataAnalysis\Paper3\';
%First Load the Flood Experiment
bd1 = load([dpath 'BottomTrack\08-03-15\HTA2_bdtrace.mat']);
fn = fieldnames(bd1);
bfile1 = dir([dpath '\BedStress\08-03-15\','*.mat']);
bfile1 = {bfile1.name};
for i = 1:3
    load([dpath '\BedStress\08-03-15\' bfile1{i}]);
    bstress1.(fn{i}) = dat;
    clear dat
end
tstart1 = datenum(2015,03,08,14,15,00);
tstop1 = datenum(2015,03,08,19,00,00);
tk1 = load([dpath 'Turbulence\08-03-15\HTA_08_TKE.mat']);
%Then load the Ebb Experiment
bd2 = load([dpath 'BottomTrack\09-03-15\HTA3_bdtrace.mat']);
bfile2 = dir([dpath '\BedStress\09-03-15\','*.mat']);
bfile2 = {bfile2.name};
for i = 1:3
    load([dpath '\BedStress\09-03-15\' bfile2{i}]);
    bstress2.(fn{i}) = dat;
    clear dat
end
tstart2 = datenum(2015,03,09,04,48,47);
tstop2 = datenum(2015,03,09,07,09,00);
tk2 = load([dpath 'Turbulence\09-03-15\HTA_09_TKE.mat']);
%Ancillary information
d50 = 1.24E-4;
psand = 0.83;pmud = 0.17;
vph = [0.242 0.271 0.240];
%Initialize Figure
f1 = figure(1);
set(f1,'PaperOrientation','portrait',...
    'position',[500 0   1000   800],...
    'renderer','painters');hold on
cl = [0.7 0.7 0.7;0.4 0.4 0.4;0 0 0];
symb = {'d','o','^'};
sp = zeros(3,1);
pl = [1 3 5;2 4 6];
for i = 1:length(fn)
    %Flood Tide
    %Plot bed trace
    sp(1) = subplot(3,2,pl(1,1));
    bdt = bd1.(fn{i}).time;
    id = find(bdt>=tstart1&bdt<=tstop1);
    bdt = bdt(id);
    bd = bd1.(fn{i}).bdist(id)-vph(i);
    plot(bdt,bd,...
        'color',cl(i,:),...
        'linewidth',1.5)
    hold on
    %Plot TKE dissipation
    sp(2) = subplot(3,2,pl(1,2));
    ept = tk1.(fn{i}).time;
    id = find(ept>=tstart1&ept<=tstop1);
    ept = ept(id);
    eps = nanmean((tk1.(fn{i}).z1.E(:,id)+tk1.(fn{i}).z2.E(:,id))./2);
    plot(ept,eps,...
        'color',cl(i,:),...
        'linewidth',1.5)
    hold on
    %Plot Bed Shear Stress
    sp(3) = subplot(3,2,pl(1,3));
    ttime = bstress1.(fn{i}).time;
    id = find(ttime>=tstart1&ttime<=tstop1);
    ttime = ttime(id);
    tmax = bstress1.(fn{i}).tmax(id);
    plot(ttime,tmax,...
        'color',cl(i,:),...
        'linewidth',1.5)
    hold on 
    if i == 3 %Critical shear stress calculation
        taucrit = critstresswu(psand,pmud,d50,'muddy');
        plot(bdt,ones(length(bdt),1)*taucrit,'--',...
            'linewidth',1.5,...
            'color','m')
    end 
    %%%
    %Ebb Tide
    %Plot bed trace
    sp(1) = subplot(3,2,pl(2,1));
    bdt = bd2.(fn{i}).time;
    id = find(bdt>=tstart2&bdt<=tstop2);
    bdt = bdt(id);
    bd = bd2.(fn{i}).bdist(id)-vph(i);
    plot(bdt,bd,...
        'color',cl(i,:),...
        'linewidth',1.5)
    hold on   
    %Plot TKE dissipation
    sp(2) = subplot(3,2,pl(2,2));
    ept = tk2.(fn{i}).time;
    id = find(ept>=tstart2&ept<=tstop2);
    ept = ept(id);
    eps = nanmean((tk2.(fn{i}).z1.E(:,id)+tk2.(fn{i}).z2.E(:,id))./2);

    plot(ept,eps,...
        'color',cl(i,:),...
        'linewidth',1.5)
    hold on
    %Plot Bed Shear Stress
    sp(3) = subplot(3,2,pl(2,3));
    ttime = bstress2.(fn{i}).time;
    id = find(ttime>=tstart2&ttime<=tstop2);
    ttime = ttime(id);
    tmax = bstress2.(fn{i}).tmax(id);
    plot(ttime,tmax,...
        'color',cl(i,:),...
        'linewidth',1.5)
    hold on
    if i == 3 %Critical shear stress calculation
        taucrit = critstresswu(psand,pmud,d50,'muddy');
        plot(bdt,ones(length(bdt),1)*taucrit,'--',...
            'linewidth',1.5,...
            'color','m')
    end  
end

%Try calculating xcorr for the lagged time series (estimate a single time
%lag for each t-s). 