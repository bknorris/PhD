%Compute the average phase angle of waves and IG bands from the
%wavelet-identified events.
clear, close all
datdir = 'd:\Mekong_W2015\DataAnalysis\Paper3\Wavelet\';
files = {'CombEventData_flood.mat';'CombEventData_ebb.mat'};
c = [207 176 126;
    60 166 74;
    4 76 41]./255;
symb = {'^';'s';'p'};
for i = 1:2
    load([datdir files{i}])
    disp([char(regexp(files{i},'(?<=\_)(.*?)(?=\.)','match')) ' tide'])
    flds = fieldnames(dat);
    figure(i)
    for j = 1:3
        disp(flds{j})
        th = deg2rad([dat.(flds{j}).wave.phase; dat.(flds{j}).ig.phase]);
        %             r = deg2rad(dat.(flds{j}).(band{k}).umag);
        %             [X,Y] = pol2cart(th,r);
        %             [alpha,~] = cart2pol(mean(X),mean(Y));
        %             fprintf('Mean circular angle:  %0.2f degrees\n',wrapTo360(rad2deg(alpha)-180))
        r = exp(1i*th);
        mu = angle(r);alpha = rad2deg(mu);
        eventl = [dat.(flds{j}).wave.eventl; dat.(flds{j}).ig.eventl];normel = eventl./max(eventl);
        bins = linspace(min(normel),max(normel),15);
        [b,n,q1,q2] = binmedian(normel,alpha,bins);
        errorbar(bins,b,q1,q2,symb{j},...
            'markerfacecolor',c(j,:),...
            'color',c(j,:),...
            'markeredgecolor','k'),hold on
    end
end