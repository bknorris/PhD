%Load the Averaged Vectrino files, and the Averaged ADV files. Find max of
%the absolute value of u,v,w and extract the height of that measurement in
%the profile. Then, divide that value by the corresponding u_inf, v_inf and
%w_inf from the ADV above. Save these data.

clear
datdir = 'd:\Mekong_W2015\DataAnalysis\Paper1\AveragedVelocities\Two\';
files = {'HTA_1Vave.mat';'HTA_2Vave.mat';'HTA_4Vave.mat';...
    'VTA_2vp1Vave.mat';'VTA_2vp2Vave.mat';'VTA_2vp3Vave.mat'};
dates = {'day1';'day2';'day3';'day4';'day4';'day4'};
VecAve = struct();
load('D:\Mekong_W2015\DataAnalysis\Paper2\BDtrack_Vels.mat')
vph = [0.062 0.063 0.061;0.2 0.2 0.2;0.5 0.5 0.5;0.07 0 0;NaN 0 0;NaN 0 0];
for i = 1:length(files)
    load([datdir files{i}])
    fn = fieldnames(Avgs);
    bn = fieldnames(data);
    if any(strcmp(dates{i},bn))
        bi = find(strcmp(dates{i},bn));
        bfn = fieldnames(data.(bn{bi}));
        bt = 1;
    else
        bt = 0;
    end
    if length(fn) > 1
        g = 1:3;
    else
        g = 1;
    end
    for ii = g
        if bt == 1
            bd = data.(bn{bi}).(bfn{ii}).bd;
            bh = vph(i,ii)-bd;
            bhadj = bh+0.002; %buffer zone of bed
        else
            bhadj = zeros(length(Avgs.(fn{ii}).time),1);
            bh = 1:length(bhadj);
        end
        vh = vph(i,ii)-0.04-linspace(0.001,0.03,35);
        VecAve.(dates{i}).(fn{ii}).time = Avgs.(fn{ii}).time;
        u = NaN(35,length(Avgs.(fn{ii}).time));
        v = NaN(35,length(Avgs.(fn{ii}).time));
        w = NaN(35,length(Avgs.(fn{ii}).time));
        uc = NaN(35,length(Avgs.(fn{ii}).time));
        %crop the bottom
        for j = 1:length(bh);
            ni = find(bhadj(j) >= vh,1,'first');
            if isempty(ni)
                ni = 35;
            end
            u(1:ni,j) = Avgs.(fn{ii}).x(1:ni,j);
            uc(1:ni,j) = Avgs.(fn{ii}).Uc(1:ni,j);
            v(1:ni,j) = Avgs.(fn{ii}).y(1:ni,j);
            w(1:ni,j) = (Avgs.(fn{ii}).z1(1:ni,j)+Avgs.(fn{ii}).z2(1:ni,j))./2;
        end
        VecAve.(dates{i}).(fn{ii}).u = u;
        VecAve.(dates{i}).(fn{ii}).v = v;
        VecAve.(dates{i}).(fn{ii}).w = w;
        VecAve.(dates{i}).(fn{ii}).Uc = uc;
    end
end
savedatdir = 'd:\Mekong_W2015\DataAnalysis\Paper2\AveragedVelocities\';
save([savedatdir 'AveVels_HTA_VTA'],'VecAve','-v7.3')

load('d:\Mekong_W2015\DataAnalysis\Paper2\ADV_avs\HTA_ADVave.mat')
cstart =  {'07-Mar-2015 15:01:00';'08-Mar-2015 15:06:00';...
    '10-Mar-2015 15:24:00';'14-Mar-2015 07:01:00'};
cstop = {'07-Mar-2015 17:07:00';'08-Mar-2015 18:03:00';...
    '10-Mar-2015 16:38:00';'14-Mar-2015 13:18:30'};
mydata = struct();fn = fieldnames(VecAve);
for i = 1:4
    if i < 4
        t1 = Vave.(fn{i}).time;
        start = datenum(cstart{i});stop = datenum(cstop{i});
        [~,id(1)] = min(abs(t1-start));
        [~,id(2)] = min(abs(t1-stop));
        uinf = Vave.(fn{i}).u(id(1):id(2));
        vinf = Vave.(fn{i}).v(id(1):id(2));
        winf = Vave.(fn{i}).w(id(1):id(2));
        dfn = fieldnames(VecAve.(fn{i}));
        for ii = 1:3
            t2 = VecAve.(fn{i}).(dfn{ii}).time;
            [~,id(1)] = min(abs(t2-start));
            [~,id(2)] = min(abs(t2-stop));
            u = VecAve.(fn{i}).(dfn{ii}).u(:,id(1):id(2));
            v = VecAve.(fn{i}).(dfn{ii}).v(:,id(1):id(2));
            w = VecAve.(fn{i}).(dfn{ii}).w(:,id(1):id(2));
            t2 = t2(id(1):id(2));
            [m,n] = size(u);
            umax = zeros(n,1);uid = zeros(n,1);
            vmax = zeros(n,1);vid = zeros(n,1);
            wmax = zeros(n,1);wid = zeros(n,1);
            for j = 1:n
                [~,ux] = max(abs(u(:,j)));
                umax(j) = u(ux,j);uid(j) = ux;
                [~,vx] = max(abs(v(:,j)));
                vmax(j) = v(vx,j);vid(j) = vx;
                [~,wx] = max(abs(w(:,j)));
                wmax(j) = w(wx,j);wid(j) = wx;
            end
            %now normalize the velocities by the over-canopy ADV vels.
            %Quality control by removing excessively large values.
            %U-vels
            uratio = umax./uinf';uratio(isnan(uratio)) = [];
            mu = mean(uratio);
            stdv = std(uratio);
            nid = find(uratio < mu+1*stdv & uratio > mu-1*stdv);
            uratio = uratio(nid);uid = uid(nid);
            utime = t2(nid);
            %V-vels
            vratio = vmax./vinf';vratio(isnan(vratio)) = [];
            mu = mean(vratio);
            stdv = std(vratio);
            nid = find(vratio < mu+1*stdv & vratio > mu-1*stdv);
            vratio = vratio(nid);vid = vid(nid);
            vtime = t2(nid);
            %W-vels
            wratio = wmax./winf';wratio(isnan(wratio)) = [];
            mu = mean(wratio);
            stdv = std(wratio);
            nid = find(wratio < mu+1*stdv & wratio > mu-1*stdv);
            wratio = wratio(nid);wid = wid(nid);
            wtime = t2(nid);
            %%%
            %Save variables to structure
            mydata.(fn{i}).(dfn{ii}).utime = utime;
            mydata.(fn{i}).(dfn{ii}).vtime = vtime;
            mydata.(fn{i}).(dfn{ii}).wtime = wtime;
            mydata.(fn{i}).(dfn{ii}).uratio = uratio;
            mydata.(fn{i}).(dfn{ii}).vratio = vratio;
            mydata.(fn{i}).(dfn{ii}).wratio = wratio;
            mydata.(fn{i}).(dfn{ii}).uid = uid;
            mydata.(fn{i}).(dfn{ii}).vid = vid;
            mydata.(fn{i}).(dfn{ii}).wid = wid;
        end
    else
        t1 = VecAve.(fn{i}).vpro3.time;
        start = datenum(cstart{i});stop = datenum(cstop{i});
        [~,id(1)] = min(abs(t1-start));
        [~,id(2)] = min(abs(t1-stop));
        uinf = nanmean(VecAve.(fn{i}).vpro3.u(:,id(1):id(2)));
        vinf = nanmean(VecAve.(fn{i}).vpro3.v(:,id(1):id(2)));
        winf = nanmean(VecAve.(fn{i}).vpro3.w(:,id(1):id(2)));
        dfn = fieldnames(VecAve.(fn{i}));
        for ii = 1:2
            t2 = VecAve.(fn{i}).(dfn{ii}).time;
            [~,id(1)] = min(abs(t2-start));
            [~,id(2)] = min(abs(t2-stop));
            u = VecAve.(fn{i}).(dfn{ii}).u(:,id(1):id(2));
            v = VecAve.(fn{i}).(dfn{ii}).v(:,id(1):id(2));
            w = VecAve.(fn{i}).(dfn{ii}).w(:,id(1):id(2));
            t2 = t2(id(1):id(2));
            [m,n] = size(u);
            umax = zeros(n,1);uid = zeros(n,1);
            vmax = zeros(n,1);vid = zeros(n,1);
            wmax = zeros(n,1);wid = zeros(n,1);
            for j = 1:n
                [~,ux] = max(abs(u(:,j)));
                umax(j) = u(ux,j);uid(j) = ux;
                [~,vx] = max(abs(v(:,j)));
                vmax(j) = v(vx,j);vid(j) = vx;
                [~,wx] = max(abs(w(:,j)));
                wmax(j) = w(wx,j);wid(j) = wx;
            end
            %now normalize the velocities by the over-canopy ADV vels.
            %Quality control by removing excessively large values.
            %U-vels
            uratio = umax./uinf';uratio(isnan(uratio)) = [];
            mu = mean(uratio);
            stdv = std(uratio);
            nid = find(uratio < mu+1*stdv & uratio > mu-1*stdv);
            uratio = uratio(nid);uid = uid(nid);
            utime = t2(nid);
            %V-vels
            vratio = vmax./vinf';vratio(isnan(vratio)) = [];
            mu = mean(vratio);
            stdv = std(vratio);
            nid = find(vratio < mu+1*stdv & vratio > mu-1*stdv);
            vratio = vratio(nid);vid = vid(nid);
            vtime = t2(nid);
            %W-vels
            wratio = wmax./winf';wratio(isnan(wratio)) = [];
            mu = mean(wratio);
            stdv = std(wratio);
            nid = find(wratio < mu+1*stdv & wratio > mu-1*stdv);
            wratio = wratio(nid);wid = wid(nid);
            wtime = t2(nid);
            %%%
            %Save variables to structure
            mydata.(fn{i}).(dfn{ii}).utime = utime;
            mydata.(fn{i}).(dfn{ii}).vtime = vtime;
            mydata.(fn{i}).(dfn{ii}).wtime = wtime;
            mydata.(fn{i}).(dfn{ii}).uratio = uratio;
            mydata.(fn{i}).(dfn{ii}).vratio = vratio;
            mydata.(fn{i}).(dfn{ii}).wratio = wratio;
            mydata.(fn{i}).(dfn{ii}).uid = uid;
            mydata.(fn{i}).(dfn{ii}).vid = vid;
            mydata.(fn{i}).(dfn{ii}).wid = wid;
        end
    end
end
%now average into 20-minute blocks
t = 20;
for i = 1:4
    dfn = fieldnames(mydata.(fn{i}));
    if i < 4
        for ii = 1:3
            t1 = 1:t:length(mydata.(fn{i}).(dfn{ii}).utime);
            urav = zeros(length(t1)-1,1);uiav = zeros(length(t1)-1,1);
            for j = 1:length(t1)-1
                urav(j) = mean(mydata.(fn{i}).(dfn{ii}).uratio(t1(j):t1(j+1)));
                uiav(j) = round(mean(mydata.(fn{i}).(dfn{ii}).uid(t1(j):t1(j+1))));
            end
            t1 = 1:t:length(mydata.(fn{i}).(dfn{ii}).vtime);
            vrav = zeros(length(t1)-1,1);viav = zeros(length(t1)-1,1);
            for j = 1:length(t1)-1
                vrav(j) = mean(mydata.(fn{i}).(dfn{ii}).vratio(t1(j):t1(j+1)));
                viav(j) = round(mean(mydata.(fn{i}).(dfn{ii}).vid(t1(j):t1(j+1))));
            end
            t1 = 1:t:length(mydata.(fn{i}).(dfn{ii}).wtime);
            wrav = zeros(length(t1)-1,1);wiav = zeros(length(t1)-1,1);
            for j = 1:length(t1)-1
                wrav(j) = mean(mydata.(fn{i}).(dfn{ii}).wratio(t1(j):t1(j+1)));
                wiav(j) = round(mean(mydata.(fn{i}).(dfn{ii}).wid(t1(j):t1(j+1))));
            end
            mydata.(fn{i}).(dfn{ii}).urav = urav;
            mydata.(fn{i}).(dfn{ii}).vrav = vrav;
            mydata.(fn{i}).(dfn{ii}).wrav = wrav;
            mydata.(fn{i}).(dfn{ii}).uiav = uiav;
            mydata.(fn{i}).(dfn{ii}).viav = viav;
            mydata.(fn{i}).(dfn{ii}).wiav = wiav;
        end
    else
        t1 = 1:t:length(mydata.(fn{i}).utime);
        urav = zeros(length(t1)-1,1);uiav = zeros(length(t1)-1,1);
        for j = 1:length(t1)-1
            urav(j) = mean(mydata.(fn{i}).uratio(t1(j):t1(j+1)));
            uiav(j) = round(mean(mydata.(fn{i}).uid(t1(j):t1(j+1))));
        end
        t1 = 1:t:length(mydata.(fn{i}).vtime);
        vrav = zeros(length(t1)-1,1);viav = zeros(length(t1)-1,1);
        for j = 1:length(t1)-1
            vrav(j) = mean(mydata.(fn{i}).vratio(t1(j):t1(j+1)));
            viav(j) = round(mean(mydata.(fn{i}).vid(t1(j):t1(j+1))));
        end
        t1 = 1:t:length(mydata.(fn{i}).wtime);
        wrav = zeros(length(t1)-1,1);wiav = zeros(length(t1),1);
        for j = 1:length(t1)-1
            wrav(j) = mean(mydata.(fn{i}).wratio(t1(j):t1(j+1)));
            wiav(j) = round(mean(mydata.(fn{i}).wid(t1(j):t1(j+1))));
        end
        mydata.(fn{i}).urav = urav;
        mydata.(fn{i}).vrav = vrav;
        mydata.(fn{i}).wrav = wrav;
        mydata.(fn{i}).uiav = uiav;
        mydata.(fn{i}).viav = viav;
        mydata.(fn{i}).wiav = wiav;
    end
end

savedatdir = 'd:\Mekong_W2015\DataAnalysis\Paper2\AveragedVelocities\';
save([savedatdir 'NormalizedVels'],'mydata','-v7.3')