clear
load('d:\Mekong_W2015\DataAnalysis\Paper3\Wavelet\CombEventData_flood.mat');
% load('g:\Data\Paper3\2015\Wavelet\CombEventData_ebb.mat')
flds = fieldnames(dat);
vars = {'eventl';'H';'usqd';'eps';'taub';'bdmed'};
for i = 1:length(flds)
    disp(['Running MLR on ' flds{i}])
    eventl = [dat.(flds{i}).wave.eventl; dat.(flds{i}).ig.eventl];
    deltbd = [dat.(flds{i}).wave.deltbd; dat.(flds{i}).ig.deltbd];
    Hs = [dat.(flds{i}).wave.sigh; dat.(flds{i}).ig.sigh];
    H = [dat.(flds{i}).wave.depth; dat.(flds{i}).ig.depth];
    uorb = [dat.(flds{i}).wave.orbwv; dat.(flds{i}).ig.orbwv];
    usqd = [dat.(flds{i}).wave.usqd; dat.(flds{i}).ig.usqd];
    taub = [dat.(flds{i}).wave.taub; dat.(flds{i}).ig.taub];
    eps = [dat.(flds{i}).wave.eps; dat.(flds{i}).ig.eps];
    id = find(deltbd>0.0015);
    data = [deltbd(id) eventl(id) H(id) usqd(id) eps(id) taub(id)];
    disp('Accretionary events')
    lma = stepwiselm(data(:,2:end),data(:,1),'constant','Upper','linear',...
        'penter',0.05,'premove',0.1,'varnames',vars)
    %%%
    id = find(deltbd<-0.0015);
    data = [deltbd(id) eventl(id) H(id) usqd(id) eps(id) taub(id)];
    disp('Erosional events')
    lme = stepwiselm(data(:,2:end),data(:,1),'constant','Upper','linear',...
        'penter',0.05,'premove',0.1,'varnames',vars)
    if i < 3
        disp('Press any key to continue')
        pause
    end
end


