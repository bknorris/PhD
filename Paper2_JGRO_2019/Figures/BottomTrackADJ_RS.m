%Adjust bottom tracking for HTA 1 and VTA (VP1): TKE FILES ONLY.

load('D:\Projects\Mekong_W2015\DataAnalysis\Paper2\BDtrack_RS.mat')
load('d:\Projects\Mekong_W2015\DataAnalysis\Paper2\RS_estimates\RStress_10min.mat')
date = {'day1';'day4'};
vph = [0.062 0.063 0.061;0.07 0 0];
for k = 1:2
    fn = fieldnames(data.(date{k}));
    dn = fieldnames(RS.(date{k}));
    if k == 2
        g = 1;
    else
        g = 1:3;
    end
    for i = g
        bd = data.(date{k}).(fn{i}).bd;
        bh = vph(k,i)-bd;
        bhadj = bh+0.002; %buffer zone of bed
        vh = vph(k,i)-0.04-linspace(0.001,0.03,30);
        for j = 1:length(bh);
            id = find(bhadj(j) >= vh,1,'first');
            RS.(date{k}).(dn{i}).uw(j,id:end) = NaN;
            RS.(date{k}).(dn{i}).vw(j,id:end) = NaN;
            RS.(date{k}).(dn{i}).uwstar(j,id:end) = NaN;
            RS.(date{k}).(dn{i}).vwstar(j,id:end) = NaN;
            RS.(date{k}).(dn{i}).ku0(j,id:end) = NaN;
            RS.(date{k}).(dn{i}).kv0(j,id:end) = NaN;
        end
    end
end
dirname = 'd:\Projects\Mekong_W2015\DataAnalysis\Paper2\RS_estimates\';
save([dirname 'RStress_10min_adapbdadj'],'RS','-v7.3')