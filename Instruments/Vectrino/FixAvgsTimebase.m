clear
dr1 = 'd:\Projects\Mekong_W2015\DataAnalysis\Paper2\AveragedVelocities\';
dr2 = 'd:\Projects\Mekong_W2015\DataAnalysis\Paper1\TKE\';
files2 = dir(dr2);files2 = {files2.name};
files1 = dir(dr1);files1 = {files1.name};
for i = 3:length(files2)
    load([dr2 files2{i}])
    load([dr1 files1{i}])
    fn = fieldnames(Stat);
    [tp,~] = size(fn);
    if tp == 1
        g = 1;
    else
        g = 1:3;
    end
    for k = g
        disp(['TKE time: ' datestr(Stat.(fn{k}).time(1))])
        disp(['Avgs time: ' datestr(Avgs.(fn{k}).time(1))])
        dt = datenum(datevec(abs(Stat.(fn{k}).time(1)-Avgs.(fn{k}).time(1))));
        if Avgs.(fn{k}).time(1) > Stat.(fn{k}).time(1)
            Avgs.(fn{k}).time = Avgs.(fn{k}).time-dt;
        elseif Avgs.(fn{k}).time(1) < Stat.(fn{k}).time(1)
            Avgs.(fn{k}).time = Avgs.(fn{k}).time+dt;
        else
            continue
        end
        disp(['New Avgs time: ' datestr(Avgs.(fn{k}).time(1))])
    end
    save([dr1 files1{i}],'Avgs','-v7.3')
    disp(['File saved: ' files1{i}])
end

