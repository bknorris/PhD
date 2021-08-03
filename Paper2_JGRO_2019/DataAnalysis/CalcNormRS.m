%Load the Averaged Vectrino files, and the Averaged ADV files. Find max of
%the absolute value of u,v,w and extract the height of that measurement in
%the profile. Then, divide that value by the corresponding u_inf, v_inf and
%w_inf from the ADV above. Save these data.

clear
load('d:\Mekong_W2015\DataAnalysis\Paper2\AveragedVelocities\AveVels_HTA_VTA.mat')
load('D:\Mekong_W2015\DataAnalysis\Paper2\RS_estimates\Final\RStress_10min_bdadj_f.mat')
load('d:\Mekong_W2015\DataAnalysis\Paper2\ADV_avs\HTA_ADVave.mat')
cstart =  {'07-Mar-2015 15:01:00';'08-Mar-2015 15:06:00';...
    '10-Mar-2015 15:24:00';'14-Mar-2015 07:01:00'};
cstop = {'07-Mar-2015 17:07:00';'08-Mar-2015 18:03:00';...
    '10-Mar-2015 16:38:00';'14-Mar-2015 13:18:30'};
mydata = struct();fn = fieldnames(RS);
for i = 1:4
    if i < 4
        t1 = Vave.(fn{i}).time;
        start = datenum(cstart{i});stop = datenum(cstop{i});
        [~,id(1)] = min(abs(t1-start));
        [~,id(2)] = min(abs(t1-stop));
        uinf = Vave.(fn{i}).u(id(1):id(2)); %just looking at x-shore
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
            uw = RS.(fn{i}).(dfn{ii}).uw(id(1):id(2),:);
            [m,n] = size(uw);
            uwmax = zeros(m,1);uid = zeros(m,1);
            for j = 1:m
                [~,ux] = max(abs(uw(j,:)));
                uwmax(j) = uw(j,ux);
                uid(j) = ux;
            end
            ustar = sqrt(abs(uwmax));
            %now average velocities to same times as ustar
            duav = zeros(length(t3)-1,1);
            for j = 1:length(t3)-1
                tid = find(t2>=t3(j) & t2<=t3(j+1));
                duav(j) = mean(abs(dU(tid)));
            end
            usnorm = ustar(1:end-1)./duav;
            %Save variables to structure
            mydata.(fn{i}).(dfn{ii}).time = t3(1:end-1);
            mydata.(fn{i}).(dfn{ii}).uid = uid(1:end-1);
            mydata.(fn{i}).(dfn{ii}).usnorm = usnorm;
            mydata.(fn{i}).(dfn{ii}).duav = duav;
            mydata.(fn{i}).(dfn{ii}).ustar = ustar(1:end-1);
        end
    else
        t1 = VecAve.day4.vpro1.time;
        start = datenum(cstart{i});stop = datenum(cstop{i});
        [~,id(1)] = min(abs(t1-start));
        [~,id(2)] = min(abs(t1-stop));
        uinf = nanmean(VecAve.day4.vpro3.u(:,id(1):id(2)));
        t2 = VecAve.day4.vpro1.time;
        [~,id(1)] = min(abs(t2-start));
        [~,id(2)] = min(abs(t2-stop));
        u = VecAve.day4.vpro1.Uc(:,id(1):id(2));
        [~,n] = size(u);
        uc = zeros(n,1);
        for j = 1:n
            [~,ux] = max(abs(u(:,j)));
            uc(j) = u(ux,j);
        end
        dU = uinf';%-uc;
        %find RS
        t3 = RS.day4.vpro1.time;
        [~,id(1)] = min(abs(t3-start));
        [~,id(2)] = min(abs(t3-stop));
        t3 = t3(id(1):id(2));
        uw = RS.day4.vpro1.uw(id(1):id(2),:);
        [m,n] = size(uw);
        uwmax = zeros(m,1);uid = zeros(m,1);
        for j = 1:m
            [~,ux] = max(abs(uw(j,:)));
            uwmax(j) = uw(j,ux);
            uid(j) = ux;
        end
        ustar = sqrt(abs(uwmax));
        %now average velocities to same times as ustar
        duav = zeros(length(t3)-1,1);
        for j = 1:length(t3)-1
            tid = find(t2>=t3(j) & t2<=t3(j+1));
            duav(j) = mean(abs(dU(tid)));
        end
        usnorm = ustar(1:end-1)./duav;
        %Save variables to structure
        mydata.(fn{i}).time = t3(1:end-1);
        mydata.(fn{i}).uid = uid(1:end-1);
        mydata.(fn{i}).usnorm = usnorm;
        mydata.(fn{i}).duav = duav;
        mydata.(fn{i}).ustar = ustar(1:end-1);
    end
end
% %now average into 20-minute blocks
% t = 20;
% for i = 1:4
%     dfn = fieldnames(mydata.(fn{i}));
%     if i < 4
%         for ii = 1:3
%             t1 = 1:t:length(mydata.(fn{i}).(dfn{ii}).utime);
%             urav = zeros(length(t1)-1,1);uiav = zeros(length(t1)-1,1);
%             for j = 1:length(t1)-1
%                 urav(j) = mean(mydata.(fn{i}).(dfn{ii}).uratio(t1(j):t1(j+1)));
%                 uiav(j) = round(mean(mydata.(fn{i}).(dfn{ii}).uid(t1(j):t1(j+1))));
%             end
%             t1 = 1:t:length(mydata.(fn{i}).(dfn{ii}).vtime);
%             vrav = zeros(length(t1)-1,1);viav = zeros(length(t1)-1,1);
%             for j = 1:length(t1)-1
%                 vrav(j) = mean(mydata.(fn{i}).(dfn{ii}).vratio(t1(j):t1(j+1)));
%                 viav(j) = round(mean(mydata.(fn{i}).(dfn{ii}).vid(t1(j):t1(j+1))));
%             end
%             t1 = 1:t:length(mydata.(fn{i}).(dfn{ii}).wtime);
%             wrav = zeros(length(t1)-1,1);wiav = zeros(length(t1)-1,1);
%             for j = 1:length(t1)-1
%                 wrav(j) = mean(mydata.(fn{i}).(dfn{ii}).wratio(t1(j):t1(j+1)));
%                 wiav(j) = round(mean(mydata.(fn{i}).(dfn{ii}).wid(t1(j):t1(j+1))));
%             end
%             mydata.(fn{i}).(dfn{ii}).urav = urav;
%             mydata.(fn{i}).(dfn{ii}).vrav = vrav;
%             mydata.(fn{i}).(dfn{ii}).wrav = wrav;
%             mydata.(fn{i}).(dfn{ii}).uiav = uiav;
%             mydata.(fn{i}).(dfn{ii}).viav = viav;
%             mydata.(fn{i}).(dfn{ii}).wiav = wiav;
%         end
%     else
%         t1 = 1:t:length(mydata.(fn{i}).utime);
%         urav = zeros(length(t1)-1,1);uiav = zeros(length(t1)-1,1);
%         for j = 1:length(t1)-1
%             urav(j) = mean(mydata.(fn{i}).uratio(t1(j):t1(j+1)));
%             uiav(j) = round(mean(mydata.(fn{i}).uid(t1(j):t1(j+1))));
%         end
%         t1 = 1:t:length(mydata.(fn{i}).vtime);
%         vrav = zeros(length(t1)-1,1);viav = zeros(length(t1)-1,1);
%         for j = 1:length(t1)-1
%             vrav(j) = mean(mydata.(fn{i}).vratio(t1(j):t1(j+1)));
%             viav(j) = round(mean(mydata.(fn{i}).vid(t1(j):t1(j+1))));
%         end
%         t1 = 1:t:length(mydata.(fn{i}).wtime);
%         wrav = zeros(length(t1)-1,1);wiav = zeros(length(t1),1);
%         for j = 1:length(t1)-1
%             wrav(j) = mean(mydata.(fn{i}).wratio(t1(j):t1(j+1)));
%             wiav(j) = round(mean(mydata.(fn{i}).wid(t1(j):t1(j+1))));
%         end
%         mydata.(fn{i}).urav = urav;
%         mydata.(fn{i}).vrav = vrav;
%         mydata.(fn{i}).wrav = wrav;
%         mydata.(fn{i}).uiav = uiav;
%         mydata.(fn{i}).viav = viav;
%         mydata.(fn{i}).wiav = wiav;
%     end
% end

savedatdir = 'd:\Mekong_W2015\DataAnalysis\Paper2\RS_estimates\';
save([savedatdir 'NormFrictionU'],'mydata','-v7.3')