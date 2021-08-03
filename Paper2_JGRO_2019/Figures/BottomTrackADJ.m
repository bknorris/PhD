%Adjust bottom tracking for HTA 1 and VTA (VP1): TKE FILES ONLY.

load('D:\Projects\Mekong_W2015\DataAnalysis\Paper2\BDtrack_TKE.mat')
load('D:\Projects\Mekong_W2015\DataAnalysis\Paper2\TKE\Vertical\VTA_2TKE.mat')
% load('D:\Projects\Mekong_W2015\DataAnalysis\Paper2\TKE\HTA_1TKE.mat')


vph = [0.07 0.063 0.061];
buffer = [0.003 0.0015 0.0015];
fn = fieldnames(data.day4);
dn = fieldnames(Stat);
for i = 1
    bd = data.day4.(fn{i}).bd;
    bh = vph(i)-bd;
    bhadj = bh+buffer(i); %buffer zone of bed
    vh = vph(i)-0.04-linspace(0.001,0.03,30);
    for j = 1:length(bh);
        id = find(bhadj(j) >= vh,1,'first');
        Stat.(dn{i}).z1.E(id:end,j) = NaN;
        Stat.(dn{i}).z2.E(id:end,j) = NaN;
%         Stat.(dn{i}).beam3.E(id:end,j) = NaN;
%         Stat.(dn{i}).beam4.E(id:end,j) = NaN;
    end
end
dirname = 'd:\Projects\Mekong_W2015\DataAnalysis\Paper2\TKE\Vertical\';
save([dirname 'VTA_2TKE_bdadj'],'Stat','-v7.3')