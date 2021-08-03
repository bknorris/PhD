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
keep = find(vifdata<5);
fprintf('Variables to keep: ')
fprintf('%s ',vals{keep})
fprintf('\n')




