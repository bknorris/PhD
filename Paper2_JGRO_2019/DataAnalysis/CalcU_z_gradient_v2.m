%Calculate velocity profile gradients for the HTA and VTA, to calculate
%Eddy viscosity profiles. 

clear
datdir = 'd:\Projects\Mekong_W2015\DataAnalysis\Paper2\AveragedVelocities\';
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
        dudz = Avgs.(fn{j}).dudz;
        vh = hab(i,j)-0.04-linspace(0.001,0.03,35);
        fn2 = fieldnames(data);
        if any(strcmp(dn{i},fn2)) && i < 5
            id = find(strcmp(dn{i},fn2));
            fn2 = fieldnames(data.(fn2{id}));
            bd = data.(dn{i}).(fn2{j}).bd;
            bh = hab(i,j)-bd;
            bhadj = bh; %buffer zone of bed
        else
            bhadj = zeros(length(dudz),1);
        end
        for k = 1:length(bhadj);
            id2 = find(bhadj(k) >= vh,1,'first');
            if isempty(id2)
                id2 = 36;
            end
            dudz(id2:end,k) = NaN;
        end
        if length(bhadj) ~= length(dudz)
            id = length(dudz)-length(bhadj);
            dudz(id2:end,end-id:end) = NaN;
        end
        dudz = dudz.*-1;
        slope.(dn{i}).(fn{j}).time = Avgs.(fn{j}).time;
        slope.(dn{i}).(fn{j}).dudz = dudz;
    end
end
dirname = 'd:\Projects\Mekong_W2015\DataAnalysis\Paper2\';
save([dirname 'U_z_gradient'],'slope','-v7.3')