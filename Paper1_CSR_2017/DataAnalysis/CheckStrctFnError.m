clear
dir1 = 'd:\Projects\Mekong_F2014\DataAnalysis\Paper1\WithN\';
dir2 = 'd:\Projects\Mekong_W2015\DataAnalysis\Paper1\WithN\';
f1 = dir([dir1 '*.mat']);f1 = {f1.name};
f2 = dir([dir2 '*.mat']);f2 = {f2.name};
%2014
horiz = [1 1 0 0 1 0 1 0];
for i = 1:length(f1)
    load([dir1 f1{i}])
    disp(f1{i})
    if horiz(i) == 1
        bin = 15;
    else
        bin = 5;
    end
    fn = fieldnames(Stat);
    if length(fn) == 1
        g = 1;
    else
        g = 1:3;
    end
    for ii = g
        dn = fieldnames(Stat.(fn{ii}));
        for j = 2:4
            E = Stat.(fn{ii}).(dn{j}).E(bin,:);
            N = Stat.(fn{ii}).(dn{j}).N(bin,:);

        
        
    