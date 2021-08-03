clear, close all
load('d:\Mekong_W2015\DataAnalysis\Paper3\Wavelet\CombEventData_flood.mat');
flds = fieldnames(dat);
vars = {'H';'Hs';'usqd';'taub';'eps';'bdmed'};
stage = {'LL';'ML';'MH';'HH'};
for i = 1:length(flds)
    disp(['Running MLR on ' flds{i}])
    exp = [dat.(flds{i}).wave.exp; dat.(flds{i}).ig.exp];
    eventl = [dat.(flds{i}).wave.eventl; dat.(flds{i}).ig.eventl];
    bdmed = [dat.(flds{i}).wave.bdmed; dat.(flds{i}).ig.bdmed];
    Hs = [dat.(flds{i}).wave.sigh; dat.(flds{i}).ig.sigh];
    H = [dat.(flds{i}).wave.depth; dat.(flds{i}).ig.depth];
    uorb = [dat.(flds{i}).wave.orbwv; dat.(flds{i}).ig.orbwv];
    usqd = [dat.(flds{i}).wave.usqd; dat.(flds{i}).ig.usqd];
    taub = [dat.(flds{i}).wave.taub; dat.(flds{i}).ig.taub];
    eps = [dat.(flds{i}).wave.eps; dat.(flds{i}).ig.eps];
    id = find(bdmed<-0.003 | bdmed > 0.003);
    data = [H(id) Hs(id) usqd(id) taub(id) eps(id)];
    
    
    
    figure(i)
    if i < 3
        disp('Press any key to continue')
        pause
    end
    %         data = [bdmed eventl H Hs usqd taub eps];
    %         lma = stepwiselm(data(:,2:end),data(:,1),'constant','Upper','linear',...
    %             'penter',0.05,'premove',0.1,'varnames',vars)
    %         if i < 3
    %             disp('Press any key to continue')
    %             pause
    %         end
end

