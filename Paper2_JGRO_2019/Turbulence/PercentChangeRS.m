%estimate % removed points for RS estimates
clear
load('D:\Projects\Mekong_W2015\DataAnalysis\Paper2\RS_estimates\RStress_10min_v5_bdadj.mat')
new = RS;clear RS
load('D:\Projects\Mekong_W2015\DataAnalysis\Paper2\RS_estimates\RStress_10min_adapbdadj.mat')
orig = RS;clear RS
fn = fieldnames(new);
for i = 1:length(fn)
    dn = fieldnames(new.(fn{i})); 
    for j = 1:3
        uwn = new.(fn{i}).(dn{j}).uw;
        uwo = orig.(fn{i}).(dn{j}).uw;
        [m,n] = size(uwn);
        pctlost = (sum(sum(isnan(uwn)))/(m*n))-(sum(sum(isnan(uwo)))/(m*n));
        fprintf('%s, %s, percentage of remaining non-nan values: %0.2f%%\n',fn{i},dn{j},100-(pctlost*100))
    end
end