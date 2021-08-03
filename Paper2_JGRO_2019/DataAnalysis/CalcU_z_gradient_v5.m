%Calculate velocity profile gradients for the HTA and VTA, to calculate
%Eddy viscosity profiles. This script is different from v3 in that it
%averages U (x) first, then calculates dudz from the entire time average. 

clear
datdir = 'D:\Projects\Mekong_W2015\DataAnalysis\Paper2\AveragedVelocities\Two\';
files = {'7March2015_Vave.mat';'8March2015_Vave.mat';'10March2015_Vave.mat';...
    '14March2015a_Vave.mat'};
load('D:\Projects\Mekong_W2015\DataAnalysis\Paper2\BottomTrack\BDtrack_Vels.mat')
dn = {'day1';'day2';'day3';'day4'};
hab = [0.062 0.063 0.061;
    0.242 0.240 0.240;
    0.550 0.550 0.550;
    0.07 0.416 0.806];
slope = struct();
cc = 1;
for i = 1:4
    load([datdir files{i}])
    disp(files{i})
    fn = fieldnames(Avgs);
    for j = 1:3
        disp(fn{j})
        u = Avgs.(fn{j}).x;
        [m,n] = size(u);
        ut = Avgs.(fn{j}).time;
        vh = hab(i,j)-0.04-linspace(0.001,0.03,m);
        fn2 = fieldnames(data);
        if any(strcmp(dn{i},fn2)) && cc < 11
            id = find(strcmp(dn{i},fn2));
            fn2 = fieldnames(data.(fn2{id}));
            bd = data.(dn{i}).(fn2{j}).bd;
            bdt = data.(dn{i}).(fn2{j}).time;
            bh = hab(i,j)-bd;
            if i == 1
                wbbl = 0.006;
            elseif i == 4
                wbbl = 0.0025; 
            end
            bhadj = bh+wbbl; %buffer zone of bed
            %need to average bdh down to the same times as u
            bd = zeros(size(ut));
            for k = 1:length(ut)-1
                idx = find(bdt>=ut(k)&bdt<=ut(k+1));
                bd(k) = nanmean(bhadj(idx));
            end
            bd(end) = bhadj(end);
        else
            bd = zeros(n,1);
        end
        uu = NaN(n,m);
        for k = 1:length(bd);
            id2 = find(bd(k) >= vh,1,'first');
            if isempty(id2)
                id2 = m;
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
        figure
        plot(umean(1:idx),(1:idx)/1000,'+r','markersize',15),hold on
        plot((1:idx)/1000,pv,'k')
        dudz = pf(1)*-1; %mult by -1 because z is positive upwards
        %calculate min, mean and max canopy velocities
        minu = min(u(1:5,:));meanu = mean(u(1:5,:));maxu = max(u(1:5,:));
        slope.(dn{i}).(fn{j}).time = Avgs.(fn{j}).time;
        slope.(dn{i}).(fn{j}).dudz = dudz;
        slope.(dn{i}).(fn{j}).minu = minu;
        slope.(dn{i}).(fn{j}).meanu  = meanu;
        slope.(dn{i}).(fn{j}).maxu = maxu;
        cc = cc+1;
    end
end
dirname = 'd:\Projects\Mekong_W2015\DataAnalysis\Paper2\AveragedVelocities\Two\';
save([dirname 'U_z_gradient_v5'],'slope','-v7.3')