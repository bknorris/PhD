clear
dr1 = 'd:\Projects\Mekong_F2014\DataAnalysis\Paper1\AveragedVelocities\';
dr2 = 'd:\Projects\Mekong_F2014\DataAnalysis\Paper1\TKE\Vertical\';
files2 = dir([dr2 '*.mat']);files2 = {files2.name};
files1 = dir([dr1 '*.mat']);files1 = {files1.name};
for i = 1:length(files2)
    load([dr2 files2{i}]),disp(files2{i})
    load([dr1 files1{i}]),disp(files1{i})
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
        if Stat.(fn{k}).time(1) > Avgs.(fn{k}).time(1)
            Stat.(fn{k}).time = Stat.(fn{k}).time-dt;
            change = 1;
        elseif Stat.(fn{k}).time(1) < Avgs.(fn{k}).time(1)
            Stat.(fn{k}).time = Stat.(fn{k}).time+dt;
            change = 1;
        else
            change = 0;
        end
        disp(['New Stat time: ' datestr(Stat.(fn{k}).time(1))])
    end
    if change == 1
    save([dr2 files2{i}],'Stat','-v7.3')
%     save([dr1 files1{i}],'Avgs','-v7.3')
    disp(['File saved: ' files2{i}])
%     disp(['File saved: ' files1{i}])
    else
        change = 0;
    end
end

