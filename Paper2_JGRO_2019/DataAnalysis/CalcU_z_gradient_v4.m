%Calculate velocity profile gradients for the HTA and VTA, to calculate
%Eddy viscosity profiles. This script is different from v3 in that it
%averages U (x) first, then calculates dudz from the entire time average. 

clear
datdir = 'd:\Mekong_W2015\DataAnalysis\Paper1\AveragedVelocities\Two\';
files = {'HTA_1Vave.mat';'HTA_2Vave.mat';'HTA_4Vave.mat';...
    'VTA_2vp1Vave.mat';'VTA_2vp2Vave.mat';'VTA_2vp3Vave.mat'};
load('D:\Mekong_W2015\DataAnalysis\Paper2\BDtrack_Vels.mat')
dn = {'day1';'day2';'day3';'day4';'day4';'day4'};
hab = [0.062 0.063 0.061;
    0.242 0.240 0.240;
    0.550 0.550 0.550;
    0.07  0 0; %zeros are dummy to get the indexing right with 6 files
    0.416 0 0;
    0.806 0 0];
slope = struct();
for i = 1:6
    load([datdir files{i}])
    disp(files{i})
    fn = fieldnames(Avgs);
    if length(fn) == 3
        g = 1:3;
    else
        g = 1;
    end
    for j = g
        disp(fn{j})
        u = Avgs.(fn{j}).x;
        vh = hab(i,j)-0.04-linspace(0.001,0.03,35);
        fn2 = fieldnames(data);
        if any(strcmp(dn{i},fn2)) && i < 5
            id = find(strcmp(dn{i},fn2));
            fn2 = fieldnames(data.(fn2{id}));
            bd = data.(dn{i}).(fn2{j}).bd;
            bh = hab(i,j)-bd;
            if i == 1
                wbbl = 0.006;
            elseif i == 4
                wbbl = 0.0025; 
            end
            bhadj = bh+wbbl; %buffer zone of bed
        else
            bhadj = zeros(length(u),1);
        end
        uu = NaN(length(u),35);
        for k = 1:length(bhadj);
            id2 = find(bhadj(k) >= vh,1,'first');
            if isempty(id2)
                id2 = 35;
            end
            uu(k,1:id2) = u(1:id2,k); %crop bottom velocities
        end
        umean = nanmean(uu);
        nid = find(isnan(umean),1,'first'); %make sure not to include any NaNs
        if isempty(nid)
            idx = 35;
        else
            idx = nid-1;
        end
        pf = polyfit(umean(1:idx),(1:idx)/1000,1);
        pv = polyval(pf,(1:idx)/1000);
%         figure
%         plot(umean(1:idx),(1:idx)/1000,'+r','markersize',15),hold on
%         plot((1:idx)/1000,pv,'k')
        dudz = pf(1)*-1; %mult by -1 because z is positive upwards
        %calculate min, mean and max canopy velocities
        minu = min(u(1:5,:));meanu = mean(u(1:5,:));maxu = max(u(1:5,:));
        slope.(dn{i}).(fn{j}).time = Avgs.(fn{j}).time;
        slope.(dn{i}).(fn{j}).dudz = dudz;
        slope.(dn{i}).(fn{j}).minu = minu;
        slope.(dn{i}).(fn{j}).meanu  = meanu;
        slope.(dn{i}).(fn{j}).maxu = maxu;
    end
end
dirname = 'd:\Mekong_W2015\DataAnalysis\Paper2\';
save([dirname 'U_z_gradient_v4'],'slope','-v7.3')