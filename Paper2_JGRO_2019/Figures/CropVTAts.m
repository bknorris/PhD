load('D:\Projects\Mekong_W2015\DataAnalysis\Paper2\TKE\Vertical\VTA_2TKE.mat')

dn = fieldnames(Stat);
for i = 1:3
        fn = fieldnames(Stat.(dn{i}));
        time = Stat.(dn{i}).time;
        start = datenum(2015,03,14,07,10,00);
        stop = datenum(2015,03,14,09,20,10);
        id = find(time >= start & time <= stop);
        Stat.(dn{i}).time = time(id);
        
        for ii = 2:3
            gn = fieldnames(Stat.(dn{i}).(fn{ii}));
            for j = 1:7
                Stat.(dn{i}).(fn{ii}).(gn{j}) = Stat.(dn{i}).(fn{ii}).(gn{j})(:,id);
            end
            for j = 8:9
                Stat.(dn{i}).(fn{ii}) = rmfield(Stat.(dn{i}).(fn{ii}),gn{j});
            end
                
        end
end
save('D:\Projects\Mekong_W2015\DataAnalysis\Paper2\TKE\Vertical\VTA_2TKE.mat','Stat','-v7.3')