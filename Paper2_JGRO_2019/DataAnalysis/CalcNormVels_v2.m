%Load the Averaged Vectrino files, and the Averaged ADV files. Find max of
%the absolute value of u,v,w and extract the height of that measurement in
%the profile. Then, divide that value by the corresponding u_inf, v_inf and
%w_inf from the ADV above. Save these data.

clear
datdir = 'd:\Projects\Mekong_W2015\DataAnalysis\Paper2\AveragedVelocities\SingleRotation\';
files = {'7March2015_Vave_sngrot.mat';'8March2015_Vave_sngrot.mat';'10March2015_Vave_sngrot.mat';...
    '14March2015a_bdadj_srot.mat'};
dates = {'day1';'day2';'day3';'day4'};
VecAve = struct();
load('d:\Projects\Mekong_W2015\DataAnalysis\Paper2\BottomTrack\BDtrack_Vels.mat')
vph = [0.062 0.063 0.061;0.2 0.2 0.2;0.5 0.5 0.5;0.07 0 0];
bt = [1 1 1;0 0 0;0 0 0;1 0 0];
for i = 1:length(files)
    load([datdir files{i}])
    fn = fieldnames(Avgs);
    bn = fieldnames(data);
    if any(strcmp(dates{i},bn))
        bi = find(strcmp(dates{i},bn));
        bfn = fieldnames(data.(bn{bi}));
    end
    for ii = 1:3
        vpt = Avgs.(fn{ii}).time;
        if bt(i,ii) == 1
            bd = data.(bn{bi}).(bfn{ii}).bd;
            bh = vph(i,ii)-bd;
            bhadj = bh+0.004; %buffer zone of bed        
            bdt = data.(bn{bi}).(bfn{ii}).time;
            newbd = zeros(length(vpt),1);
            for j = 1:length(vpt)-1
            iid = find(bdt>=vpt(j)&bdt<=vpt(j+1));
            newbd(j) = nanmean(bhadj(iid));
            end
            newbd(end) = bhadj(end);
            bhadj = newbd;
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
        uw = NaN(35,length(Avgs.(fn{ii}).time));
        %crop the bottom
        for j = 1:length(bhadj);
            ni = find(bhadj(j) >= vh,1,'first');
            if isempty(ni)
                ni = 35;
            end
            u(1:ni,j) = Avgs.(fn{ii}).x(1:ni,j);
            uc(1:ni,j) = Avgs.(fn{ii}).Uc(1:ni,j);
            uw(1:ni,j) = Avgs.(fn{ii}).Uw(1:ni,j);
            v(1:ni,j) = Avgs.(fn{ii}).y(1:ni,j);
            w(1:ni,j) = (Avgs.(fn{ii}).z1(1:ni,j)+Avgs.(fn{ii}).z2(1:ni,j))./2;
        end
        VecAve.(dates{i}).(fn{ii}).u = u;
        VecAve.(dates{i}).(fn{ii}).v = v;
        VecAve.(dates{i}).(fn{ii}).w = w;
        VecAve.(dates{i}).(fn{ii}).Uc = uc;
        VecAve.(dates{i}).(fn{ii}).Uw = uw;
    end
end
savedatdir = 'd:\Projects\Mekong_W2015\DataAnalysis\Paper2\AveragedVelocities\';
save([savedatdir 'AveVels_HTA_VTA'],'VecAve','-v7.3')
break
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
        for ii = 1:2
            t1 = VecAve.(fn{i}).(dfn{ii}).time;
            start = datenum(cstart{i});stop = datenum(cstop{i});
            [~,id(1)] = min(abs(t1-start));
            [~,id(2)] = min(abs(t1-stop));
            uinf = nanmean(VecAve.(fn{i}).vpro3.u(:,id(1):id(2)));
            vinf = nanmean(VecAve.(fn{i}).vpro3.v(:,id(1):id(2)));
            winf = nanmean(VecAve.(fn{i}).vpro3.w(:,id(1):id(2)));
            dfn = fieldnames(VecAve.(fn{i}));
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
%now average into 10-minute blocks
t = 10;
for i = 1:4
    dfn = fieldnames(mydata.(fn{i}));
    for ii = 1:length(dfn)
        t1 = 1:t:length(mydata.(fn{i}).(dfn{ii}).utime);
        urav = zeros(length(t1)-1,1);uiav = zeros(length(t1)-1,1);ut = zeros(length(t1)-1,1);
        for j = 1:length(t1)-1
            urav(j) = mean(mydata.(fn{i}).(dfn{ii}).uratio(t1(j):t1(j+1)));
            uiav(j) = round(mean(mydata.(fn{i}).(dfn{ii}).uid(t1(j):t1(j+1))));
            ut(j) = mydata.(fn{i}).(dfn{ii}).utime(t1(j));
        end
        t1 = 1:t:length(mydata.(fn{i}).(dfn{ii}).vtime);vt = zeros(length(t1)-1,1);
        vrav = zeros(length(t1)-1,1);viav = zeros(length(t1)-1,1);
        for j = 1:length(t1)-1
            vrav(j) = mean(mydata.(fn{i}).(dfn{ii}).vratio(t1(j):t1(j+1)));
            viav(j) = round(mean(mydata.(fn{i}).(dfn{ii}).vid(t1(j):t1(j+1))));
            vt(j) = mydata.(fn{i}).(dfn{ii}).vtime(t1(j));
        end
        t1 = 1:t:length(mydata.(fn{i}).(dfn{ii}).wtime);wt = zeros(length(t1)-1,1);
        wrav = zeros(length(t1)-1,1);wiav = zeros(length(t1)-1,1);
        for j = 1:length(t1)-1
            wrav(j) = mean(mydata.(fn{i}).(dfn{ii}).wratio(t1(j):t1(j+1)));
            wiav(j) = round(mean(mydata.(fn{i}).(dfn{ii}).wid(t1(j):t1(j+1))));
            wt(j) = mydata.(fn{i}).(dfn{ii}).wtime(t1(j));
        end
        mydata.(fn{i}).(dfn{ii}).utime2 = ut;
        mydata.(fn{i}).(dfn{ii}).vtime2 = vt;
        mydata.(fn{i}).(dfn{ii}).wtime2 = wt;
        mydata.(fn{i}).(dfn{ii}).urav = urav;
        mydata.(fn{i}).(dfn{ii}).vrav = vrav;
        mydata.(fn{i}).(dfn{ii}).wrav = wrav;
        mydata.(fn{i}).(dfn{ii}).uiav = uiav;
        mydata.(fn{i}).(dfn{ii}).viav = viav;
        mydata.(fn{i}).(dfn{ii}).wiav = wiav;
    end
end

savedatdir = 'd:\Mekong_W2015\DataAnalysis\Paper2\AveragedVelocities\';
save([savedatdir 'NormalizedVels'],'mydata','-v7.3')