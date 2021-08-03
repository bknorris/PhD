%Load the Averaged Vectrino files, and the Averaged ADV files. Find max of
%the absolute value of u,v,w and extract the height of that measurement in
%the profile. Then, divide that value by the corresponding u_inf, v_inf and
%w_inf from the ADV above. Save these data.

clear
load('d:\Projects\Mekong_W2015\DataAnalysis\Paper2\AveragedVelocities\AveVels_HTA_VTA.mat')
load('d:\Projects\Mekong_W2015\DataAnalysis\Paper2\RS_estimates\RStress_10min_v5_bdadj.mat')
load('d:\Projects\Mekong_W2015\DataAnalysis\Paper2\AveragedVelocities\HTA_ADVave.mat')
cstart =  {'07-Mar-2015 15:01:00';'08-Mar-2015 15:06:00';...
    '10-Mar-2015 15:24:00';'14-Mar-2015 04:40:00'};
cstop = {'07-Mar-2015 17:07:00';'08-Mar-2015 18:03:00';...
    '10-Mar-2015 16:38:00';'14-Mar-2015 10:40:30'};
mydata = struct();fn = fieldnames(RS);
for i = 1:4
    if i < 4
        t1 = Vave.(fn{i}).time;
        start = datenum(cstart{i});stop = datenum(cstop{i});
        [~,id(1)] = min(abs(t1-start));
        [~,id(2)] = min(abs(t1-stop));
        uinf = Vave.(fn{i}).Uw(id(1):id(2)); %just looking at x-shore
        dfn = fieldnames(VecAve.(fn{i}));
        for ii = 1:3
            t2 = VecAve.(fn{i}).(dfn{ii}).time;
            [~,id(1)] = min(abs(t2-start));
            [~,id(2)] = min(abs(t2-stop));
            t2 = t2(id(1):id(2));
            u = VecAve.(fn{i}).(dfn{ii}).Uc(:,id(1):id(2));
            [~,n] = size(u);
            uc = zeros(n,1);
            for j = 1:n
                [~,ux] = max(abs(u(:,j)));
                uc(j) = u(ux,j);
            end
            dU = uinf';%-uc;
            %find RS
            t3 = RS.(fn{i}).(dfn{ii}).time;
            [~,id(1)] = min(abs(t3-start));
            [~,id(2)] = min(abs(t3-stop));
            t3 = t3(id(1):id(2));
            uw = nanmean(RS.(fn{i}).(dfn{ii}).uw(id(1):id(2),:),2);
            vw = nanmean(RS.(fn{i}).(dfn{ii}).uw(id(1):id(2),:),2);
            ustar = sqrt(uw.^2+vw.^2);
            %now average velocities to same times as ustar
            duav = zeros(length(t3)-1,1);
            for j = 1:length(t3)-1
                tid = find(t2>=t3(j) & t2<=t3(j+1));
                duav(j) = mean(abs(dU(tid)));
            end
            usnorm = (ustar(1:end-1))./duav;
            %Save variables to structure
            mydata.(fn{i}).(dfn{ii}).time = t3(1:end-1);
            mydata.(fn{i}).(dfn{ii}).usnorm = usnorm;
            mydata.(fn{i}).(dfn{ii}).duav = duav;
            mydata.(fn{i}).(dfn{ii}).ustar = ustar(1:end-1);
        end
    else
        for ii = 1:3
            dfn = fieldnames(VecAve.(fn{i}));
            t1 = VecAve.(fn{i}).(dfn{ii}).time;
            start = datenum(cstart{i});stop = datenum(cstop{i});
            [~,id(1)] = min(abs(t1-start));
            [~,id(2)] = min(abs(t1-stop));
            uinf = nanmean(VecAve.(fn{i}).vpro3.Uw(:,id(1):id(2))); %just looking at x-shore
            t2 = VecAve.(fn{i}).(dfn{ii}).time;
            [~,id(1)] = min(abs(t2-start));
            [~,id(2)] = min(abs(t2-stop));
            t2 = t2(id(1):id(2));
            u = VecAve.(fn{i}).(dfn{ii}).Uc(:,id(1):id(2));
            [~,n] = size(u);
            uc = zeros(n,1);
            for j = 1:n
                [~,ux] = max(abs(u(:,j)));
                uc(j) = u(ux,j);
            end
            dU = uinf';%-uc;
            %find RS
            t3 = RS.(fn{i}).(dfn{ii}).time;
            [~,id(1)] = min(abs(t3-start));
            [~,id(2)] = min(abs(t3-stop));
            t3 = t3(id(1):id(2));
            uw = nanmean(RS.(fn{i}).(dfn{ii}).uw(id(1):id(2),:),2);
            vw = nanmean(RS.(fn{i}).(dfn{ii}).uw(id(1):id(2),:),2);
            ustar = sqrt(uw.^2+vw.^2);
            %now average velocities to same times as ustar
            duav = zeros(length(t3)-1,1);
            for j = 1:length(t3)-1
                tid = find(t2>=t3(j) & t2<=t3(j+1));
                duav(j) = mean(abs(dU(tid)));
            end
            usnorm = (ustar(1:end-1))./duav;
            %Save variables to structure
            mydata.(fn{i}).(dfn{ii}).time = t3(1:end-1);
            mydata.(fn{i}).(dfn{ii}).usnorm = usnorm;
            mydata.(fn{i}).(dfn{ii}).duav = duav;
            mydata.(fn{i}).(dfn{ii}).ustar = ustar(1:end-1);
        end
    end
end
for i = 1:4
    dfn = fieldnames(mydata.(fn{i}));
    for ii = 1:3
        xs = nanmean(mydata.(fn{i}).(dfn{ii}).usnorm);
        disp([fn{i} ' ' dfn{ii} ' drag coef: ' num2str(xs)])
    end
end
% savedatdir = 'd:\Projects\Mekong_W2015\DataAnalysis\Paper2\RS_estimates\';
% save([savedatdir 'NormFrictionU_v2'],'mydata','-v7.3')