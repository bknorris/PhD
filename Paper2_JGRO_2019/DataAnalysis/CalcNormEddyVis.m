%Normalize the eddy viscosity by z*deltaU where deltaU is the difference
%between the over-canopy velocity, u_h_c and the in-canopy velocity, u_c.
clear
load('d:\Mekong_W2015\DataAnalysis\Paper2\RS_estimates\Final\RStress_10min_bdadj_f.mat')
load('D:\Mekong_W2015\DataAnalysis\Paper2\U_z_gradient_v4.mat')
load('d:\Mekong_W2015\DataAnalysis\Paper2\ADV_avs\HTA_ADVave.mat')
cstart =  {'07-Mar-2015 15:01:00';'08-Mar-2015 15:06:00';...
    '10-Mar-2015 15:24:00';'14-Mar-2015 07:01:00'};
cstop = {'07-Mar-2015 17:07:00';'08-Mar-2015 18:03:00';...
    '10-Mar-2015 16:38:00';'14-Mar-2015 13:18:30'};
hc = [0.64,0.64,0.64,0.6];
%first, calculate nu
fn = fieldnames(RS);
for i = 1:4
    dn = fieldnames(RS.(fn{i}));
    for ii = 1:3
        uw = RS.(fn{i}).(dn{ii}).uw;
        dudz = repmat(slope.(fn{i}).(dn{ii}).dudz,1,35);
        %eddy viscosity is Re stress/velocity gradient
        [m,n] = size(uw);nue = zeros(m,n);
        for j = 1:m
            nue(j,:) = -1*(uw(j,:)./dudz);
        end
        nue(isinf(nue)) = NaN;
        %save to file
        nu.(fn{i}).(dn{ii}).time = RS.(fn{i}).(dn{ii}).time;
        nu.(fn{i}).(dn{ii}).nu = nue;
    end
end
%next, calulate z*deltaU
for i = 1:4
    if i < 4
        t1 = Vave.(fn{i}).time;
        start = datenum(cstart{i});stop = datenum(cstop{i});
        [~,id(1)] = min(abs(t1-start));
        [~,id(2)] = min(abs(t1-stop));
        uinf = Vave.(fn{i}).u(id(1):id(2)); %just looking at x-shore
        dfn = fieldnames(slope.(fn{i}));
        for ii = 1:3
            t2 = slope.(fn{i}).(dfn{ii}).time;
            [~,id(1)] = min(abs(t2-start));
            [~,id(2)] = min(abs(t2-stop));
            t2 = t2(id(1):id(2));
            u = slope.(fn{i}).(dfn{ii}).meanu(:,id(1):id(2));
            dU = uinf-u;
            %now average velocities to same times as nu
            t3 = nu.(fn{i}).(dfn{ii}).time;
            [~,id(1)] = min(abs(t3-start));
            [~,id(2)] = min(abs(t3-stop));
            t3 = t3(id(1):id(2));
            duav = zeros(length(t3)-1,1);
            for j = 1:length(t3)-1
                tid = find(t2>=t3(j) & t2<=t3(j+1));
                duav(j) = mean(dU(tid));
            end
            %find max(nu)
            nue = nu.(fn{i}).(dfn{ii}).nu(id(1):id(2),:);
            [m,n] = size(nue);
            nmax = zeros(m,1);nid = zeros(m,1);
            for j = 1:m
                [~,ux] = max(abs(nue(j,:)));
                nmax(j) = nue(j,ux);
                nid(j) = ux;
            end
            evnorm = nmax(1:end-1)./(duav);
            %quality control
            evnorm(isnan(evnorm)) = [];
            mu = mean(evnorm);
            stdv = std(evnorm);
            nn = find(evnorm < mu+0.5*stdv & evnorm > mu-0.5*stdv);
            evnorm = evnorm(nn);nid = nid(nn);
            time = t3(nn);
            %save variables to structure
            mydata.(fn{i}).(dfn{ii}).time = time;
            mydata.(fn{i}).(dfn{ii}).nid = nid;
            mydata.(fn{i}).(dfn{ii}).nmax = nmax(nn);
            mydata.(fn{i}).(dfn{ii}).evnorm = evnorm;
            mydata.(fn{i}).(dfn{ii}).duav = duav(nn);
            mydata.(fn{i}).(dfn{ii}).hc = hc(i);
        end
    else
        for ii = 1:2            
            dfn = fieldnames(slope.(fn{i}));
            t1 = slope.(fn{i}).(dfn{ii}).time;
            start = datenum(cstart{i});stop = datenum(cstop{i});
            [~,id(1)] = min(abs(t1-start));
            [~,id(2)] = min(abs(t1-stop));
            uinf = slope.(fn{i}).vpro3.meanu(id(1):id(2)); %just looking at x-shore
            t2 = slope.(fn{i}).(dfn{ii}).time;
            [~,id(1)] = min(abs(t2-start));
            [~,id(2)] = min(abs(t2-stop));
            t2 = t2(id(1):id(2));
            u = slope.(fn{i}).(dfn{ii}).meanu(:,id(1):id(2));
            dU = uinf-u;
            %now average velocities to same times as nu
            t3 = nu.(fn{i}).(dfn{ii}).time;
            [~,id(1)] = min(abs(t3-start));
            [~,id(2)] = min(abs(t3-stop));
            t3 = t3(id(1):id(2));
            duav = zeros(length(t3)-1,1);
            for j = 1:length(t3)-1
                tid = find(t2>=t3(j) & t2<=t3(j+1));
                duav(j) = mean(dU(tid));
            end
            %find max(nu)
            nue = nu.(fn{i}).(dfn{ii}).nu(id(1):id(2),:);
            [m,n] = size(nue);
            nmax = zeros(m,1);nid = zeros(m,1);
            for j = 1:m
                [~,ux] = max(abs(nue(j,:)));
                nmax(j) = nue(j,ux);
                nid(j) = ux;
            end
            evnorm = nmax(1:end-1)./(duav);
            %quality control
            evnorm(isnan(evnorm)) = [];
            mu = mean(evnorm);
            stdv = std(evnorm);
            nn = find(evnorm < mu+0.5*stdv & evnorm > mu-0.5*stdv);
            evnorm = evnorm(nn);nid = nid(nn);
            time = t3(nn);
            %save variables to structure
            mydata.(fn{i}).(dfn{ii}).time = time;
            mydata.(fn{i}).(dfn{ii}).nid = nid;
            mydata.(fn{i}).(dfn{ii}).nmax = nmax(nn);
            mydata.(fn{i}).(dfn{ii}).evnorm = evnorm;
            mydata.(fn{i}).(dfn{ii}).duav = duav(nn);
            mydata.(fn{i}).(dfn{ii}).hc = hc(i);
        end
    end
end

savedatdir = 'd:\Mekong_W2015\DataAnalysis\Paper2\RS_estimates\';
save([savedatdir 'NormEddyVis'],'mydata','-v7.3')
