%Adjust bottom tracking for HTA 1 and VTA (VP1): AVGS files only.

load('D:\Projects\Mekong_W2015\DataAnalysis\Paper1\BDTrack\F2F2_BDtrack_TKE.mat')
load('d:\Projects\Mekong_W2015\DataAnalysis\Paper1\AveragedVelocities\F2F2_2Vave.mat')
% load('D:\Projects\Mekong_W2015\DataAnalysis\Paper2\TKE\HTA_1TKE.mat')


vph = [0.063 0.065 0.065];
buffer = [0.0025 0.0025 0.0018];
fn = fieldnames(data);
dn = fieldnames(Stat);
for i = 1:3
    bd = data.(fn{i}).bd;
    bh = vph(i)-bd;
    bhadj = bh+buffer(i); %buffer zone of bed
    vh = vph(i)-0.04-linspace(0.001,0.03,30);
    for j = 1:length(bh);
        id = find(bhadj(j) >= vh,1,'first');
        Avgs.(dn{i}).x(id:end,j) = NaN;
        Avgs.(dn{i}).y(id:end,j) = NaN;
        Avgs.(dn{i}).z1(id:end,j) = NaN;
        Avgs.(dn{i}).z2(id:end,j) = NaN;
    end
end
dirname = 'd:\Projects\Mekong_W2015\DataAnalysis\Paper1\AveragedVelocities\';
save([dirname 'F2F2_2Avgs_bdadj'],'Avgs','-v7.3')