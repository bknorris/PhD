%Adjust bottom tracking for HTA 1 and VTA (VP1): AVG FILES ONLY.
clear
load('d:\Projects\Mekong_W2015\DataAnalysis\Paper2\BottomTrack\BDtrack_Vels.mat')
files = dir(['d:\Projects\Mekong_W2015\DataAnalysis\Paper2\AveragedVelocities\SingleRotation\' '*.mat']);    
sdatdir = 'd:\Projects\Mekong_W2015\DataAnalysis\Paper2\AveragedVelocities\SingleRotation\';
files = {files.name};
ord = [4];
vph = [0.065 0.065 0.065;0.07 0 0];
date = {'day1';'day4'};
for k = 1
    dn = fieldnames(data.(date{k}));
    load(['d:\Projects\Mekong_W2015\DataAnalysis\Paper2\AveragedVelocities\SingleRotation\' files{ord(k)}])
    for i = 1:3
        if k == 2 && i > 1
            continue
        else
            bd = data.(date{k}).(dn{i}).bd;
            bdt = data.(date{k}).(dn{i}).time;
            bh = vph(k,i)-bd;
            if i == 1
                bz = 0.006;
            else
                bz = 0.002;
            end
            bhadj = bh+bz; %buffer zone of bed
            vh = vph(k,i)-0.04-linspace(0.001,0.03,35);
            vpt = Avgs.(dn{i}).time;
            bhadj2 = zeros(length(vpt),1);
            for j = 1:length(vpt)-1
                tid = find(bdt>=vpt(j)&bdt<=vpt(j+1));
                bhadj2(j) = mean(bhadj(tid));
            end
            bhadj2(end) = bh(end);
            for j = 1:length(bhadj2);
                id = find(bhadj2(j) >= vh,1,'first');
                Avgs.(dn{i}).x(id:end,j) = NaN;
                Avgs.(dn{i}).y(id:end,j) = NaN;
                Avgs.(dn{i}).z1(id:end,j) = NaN;
                Avgs.(dn{i}).z2(id:end,j) = NaN;
                Avgs.(dn{i}).Umag(id:end,j) = NaN;
                Avgs.(dn{i}).ucubed(id:end,j) = NaN;
            end
        end
    end
    vname = regexprep(files{ord(k)},'Vave_sngrot.mat','bdadj_srot');
    save([sdatdir vname],'Avgs','-v7.3')
end
