%Calculate velocity profile gradients for the HTA and VTA, to calculate
%Eddy viscosity profiles. 

clear
datdir = 'd:\Projects\Mekong_W2015\DataAnalysis\Paper1\AveragedVelocities\Two\';
files = {'HTA_1Vave.mat';'HTA_2Vave.mat';'HTA_4Vave.mat';...
    'VTA_2vp1Vave.mat';'VTA_2vp2Vave.mat';'VTA_2vp3Vave.mat'};
load('D:\Projects\Mekong_W2015\DataAnalysis\Paper2\BDtrack_Vels.mat')
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
            bhadj = bh+0.008; %buffer zone of bed
        else
            bhadj = zeros(length(u),1);
        end
        dudz = NaN(length(u),35);
        for k = 1:length(bhadj);
            id2 = find(bhadj(k) >= vh,1,'first');
            if isempty(id2)
                id2 = 36;
            end
            uu = u(:,k);
            if length(uu) < 2 || sum(uu) == 0;
                dudz(k,1:id2-1) = NaN;
            else
                ind = [1 3:2:35];
                ind2 = 2:2:35;
                grad = NaN(1,35);
                for jj = 1:length(ind)-1
                    b = gradient([uu(ind(jj)) uu(ind(jj+1))],0.003);
                    grad(ind2(jj)) = b(1);
                end
                grad = grad*-1;
%                 grad = fixgaps(grad);
                grad(id2:end) = NaN;
                dudz(k,:) = grad;
            end 
        end
        dUdZ = NaN(length(u),35);
        for k = 1:length(u)
            %calculate median absolute deviation for outlier id
            b = dudz(k,:);
            mad = nanmedian(abs(b-nanmedian(b)));
            id = find(b > 3*mad | b < -3*mad);b(id) = NaN;
            dUdZ(k,:) = b;
        end
        %calculate min, mean and max canopy velocities
        minu = min(u(1:5,:));meanu = mean(u(1:5,:));maxu = max(u(1:5,:));
        slope.(dn{i}).(fn{j}).time = Avgs.(fn{j}).time;
        slope.(dn{i}).(fn{j}).dudz = dUdZ;
        slope.(dn{i}).(fn{j}).minu = minu;
        slope.(dn{i}).(fn{j}).meanu  = meanu;
        slope.(dn{i}).(fn{j}).maxu = maxu;
    end
end
dirname = 'd:\Projects\Mekong_W2015\DataAnalysis\Paper2\';
save([dirname 'U_z_gradient'],'slope','-v7.3')