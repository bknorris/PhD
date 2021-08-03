clear,close all
load('d:\Mekong_W2015\DataAnalysis\Paper3\Wavelet\CombEventData_flood.mat');
flds = fieldnames(dat);
band = {'wave';'ig'};
vars = {'eventl';'H';'usqd';'eps';'taub';'bdmed'};
for i = 1:length(flds)
    disp(['Running MLR on ' flds{i}])
    for ii = 1:2
        disp([band{ii} ' band'])
        exp = dat.(flds{i}).(band{ii}).exp;
        eventl = dat.(flds{i}).(band{ii}).eventl;
        bdmed = dat.(flds{i}).(band{ii}).bdmed;
        Hs = dat.(flds{i}).(band{ii}).sigh;
        H = dat.(flds{i}).(band{ii}).depth;
        uorb = dat.(flds{i}).(band{ii}).orbwv;
        usqd = dat.(flds{i}).(band{ii}).usqd;
        taub = dat.(flds{i}).(band{ii}).taub;
        eps = dat.(flds{i}).(band{ii}).eps;
        data = [bdmed eventl H usqd taub eps];
        lma = stepwiselm(data(:,2:end),data(:,1),'constant','Upper','linear',...
            'penter',0.05,'premove',0.1,'varnames',vars)
        if i < 3
            disp('Press any key to continue')
            pause
        end
    end
end
