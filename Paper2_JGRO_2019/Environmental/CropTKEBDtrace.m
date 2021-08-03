%Adjust bottom tracking for HTA 1 and VTA (VP1): TKE FILES ONLY.
clear
load('d:\Projects\Mekong_W2015\DataAnalysis\Paper2\BottomTrack\BDtrack_Vels.mat')
files = dir(['d:\Projects\Mekong_W2015\DataAnalysis\Paper2\Turbulence\' '*.mat']);    
sdatdir = 'd:\Projects\Mekong_W2015\DataAnalysis\Paper2\Turbulence\';
files = {files.name};
ord = [2 5];
vph = [0.062 0.063 0.061;0.07 0 0];
date = {'day1';'day4'};
for k = 1:2
    dn = fieldnames(data.(date{k}));
    load(['d:\Projects\Mekong_W2015\DataAnalysis\Paper2\Turbulence\' files{ord(k)}])
    for i = 1:3
        if k == 2 && i > 1
            continue
        else
            bd = data.(date{k}).(dn{i}).bd;
            bdt = data.(date{k}).(dn{i}).time;
            bh = vph(k,i)-bd;
            bhadj = bh+0.004; %buffer zone of bed
            vh = vph(k,i)-0.04-linspace(0.001,0.03,35);
            vpt = Stat.(dn{i}).time;
            bhadj2 = zeros(length(vpt),1);
            for j = 1:length(vpt)-1
                tid = find(bdt>=vpt(j)&bdt<=vpt(j+1));
                bhadj2(j) = mean(bhadj(tid));
            end
            bhadj2(end) = bh(end);
            for j = 1:length(bhadj2);
                id = find(bhadj2(j) >= vh,1,'first');
                Stat.(dn{i}).z1.E(id:end,j) = NaN;
                Stat.(dn{i}).z2.E(id:end,j) = NaN;
            end
        end
    end
    vname = regexprep(files{ord(k)},'VelsTKE.mat','TKEbdadj');
    save([sdatdir vname],'Stat','-v7.3')
end
