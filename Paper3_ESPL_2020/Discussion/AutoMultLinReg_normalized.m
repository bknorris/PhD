%MLR using zscore to normalize input values
clear, close all
load('e:\Mekong_W2015\DataAnalysis\Paper3\Wavelet\CombEventData_ebb_v2.mat');
flds = fieldnames(dat);
vals = {'H';'Hs';'uorb';'umag';'taub';'eps';'bdmed'};
stage = {'LL';'ML';'MH';'HH'};c = jet(3);
figno = [3 4 5;6 7 8];
for i = 1:length(flds)
    disp(['Running MLR on ' flds{i}])
    exp = [dat.(flds{i}).wave.exp; dat.(flds{i}).ig.exp];
    eventl = [dat.(flds{i}).wave.eventl; dat.(flds{i}).ig.eventl];
    bdmed = [dat.(flds{i}).wave.bdmed; dat.(flds{i}).ig.bdmed];
    Hs = [dat.(flds{i}).wave.sigh; dat.(flds{i}).ig.sigh];
    H = [dat.(flds{i}).wave.depth; dat.(flds{i}).ig.depth];
    uorb = [dat.(flds{i}).wave.orbwv; dat.(flds{i}).ig.orbwv];
    umag = [dat.(flds{i}).wave.umag; dat.(flds{i}).ig.umag];
    taub = [dat.(flds{i}).wave.taub; dat.(flds{i}).ig.taub];
    eps = [dat.(flds{i}).wave.eps; dat.(flds{i}).ig.eps];
    
    %Check for multi-collinearity!
    vals = {'H';'Hs';'uorb';'umag';'taub';'eps';'bdmed'};
    combinations = [0 1 1 1 1 1;
        1 0 1 1 1 1;...
        1 1 0 1 1 1;...
        1 1 1 0 1 1;...
        1 1 1 1 0 1;...
        1 1 1 1 1 0];
    alldata = [H Hs uorb umag taub eps];
    vifdata = zeros(1,6);
    for j = 1:6 %combinations
        row = combinations(j,:);
        dep = find(row==0);
        indep = find(row~=0);
        [~,bint,~,~,~]  = regress(zscore(alldata(:,dep)),[ones(size(H)) zscore(alldata(:,indep))],0.01);
        %check to see which variables had the highest variance
        bint = abs(mean(bint,2));
        hiv = find(bint<0.5);
        %remove these values from the rest and re-run regress
        newrow = [row' bint<0.5];newrow = newrow(:,1).*newrow(:,2);
        [b,bint,r,rint,stats]  = regress(zscore(alldata(:,dep)),[ones(size(H)) zscore(alldata(:,newrow~=0))],0.01);
        VIF = 1/(1-stats(1));
        notincl = vals(dep);incl = vals(newrow~=0);
        fprintf('VIF of %s vs ',notincl{:})
        fprintf('%s ',incl{:})
        fprintf(': %0.3f\n',VIF)
        vifdata(j) = VIF;
    end
    %reduce number of variables by the VIF
    keep = find(vifdata<2);
    fprintf('Variables to keep: ')  
    fprintf('%s ',vals{keep})
    fprintf('\n')
    lma = stepwiselm(zscore(alldata(:,keep)),zscore(bdmed),'constant','Upper','linear',...
        'penter',0.01,'premove',0.1,'varnames',[vals(keep)' {'bdmed'}])
    pause
end


