clear
datdir = 'd:\Projects\Mekong_W2015\DataAnalysis\Paper2\AveragedVelocities\Two\';
files = {'7March2015_Vave.mat';'8March2015_Vave.mat';'10March2015_Vave.mat';...
    '14March2015a_Vave.mat'};
load('d:\Projects\Mekong_W2015\DataAnalysis\Paper2\RS_estimates\RStress_10min_v5_bdadj.mat')
cstart =  {'07-Mar-2015 15:01:00';'08-Mar-2015 15:06:00';...
    '10-Mar-2015 15:24:00';'14-Mar-2015 04:40:00'};
cstop = {'07-Mar-2015 17:07:00';'08-Mar-2015 18:03:00';...
    '10-Mar-2015 16:38:00';'14-Mar-2015 10:40:30'};
mydata = struct();fn = fieldnames(RS);
for i = 1:4
    load([datdir files{i}])
    dfn = fieldnames(Avgs);        
    for ii = 1:3
        start = datenum(cstart{i});stop = datenum(cstop{i});
        %compute mean velocities at instrument
        dU = Avgs.(dfn{ii}).Umag;
        %now find mean RS
        uw = nanmean(RS.(fn{i}).(dfn{ii}).uw,2);
        vw = nanmean(RS.(fn{i}).(dfn{ii}).vw,2);
        ustar = sqrt(uw.^2+vw.^2);
        usnorm = ustar./nanmean(dU.^2)';
        %Save variables to structure
        mydata.(fn{i}).(dfn{ii}).time = Avgs.(dfn{ii}).time;
        mydata.(fn{i}).(dfn{ii}).uid = repmat(15,length(uw),1);
        mydata.(fn{i}).(dfn{ii}).usnorm = usnorm;
        mydata.(fn{i}).(dfn{ii}).duav = dU';
        mydata.(fn{i}).(dfn{ii}).ustar = ustar(1:end-1);
    end
end
savedatdir = 'd:\Projects\Mekong_W2015\DataAnalysis\Paper2\RS_estimates\';
save([savedatdir 'NormFrictionU_v5'],'mydata','-v7.3')