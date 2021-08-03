%Calculate Uw and Uc, ST #, & Re # for the HTA and VTA.
clear
datdir = 'd:\Projects\Mekong_W2015\DataAnalysis\Paper1\';
files = {'HTA_1Vels.mat';'HTA_2Vels.mat';'HTA_4Vels.mat';...
    'VTA_2vp1Vels.mat';'VTA_2vp2Vels.mat';'VTA_2vp3Vels.mat'};
days = {'day1';'day2';'day3';'day4';'day4';'day4'};
heading = [20 20 20 96 96 96];
d = [0.019 0.019 0.019;
    0.018 0.018 0.018;
    0.008 0.008 0.008;
    0.013 NaN NaN
    0.004 NaN NaN
    0 NaN NaN];
fpeaks = [2.3 2.5 NaN;
    2.03 2.16 2.4;
    1.85 NaN NaN];
data = struct();
for i = 1:6
    disp(['Loading ' files{i}])
    load([datdir files{i}])
    fn = fieldnames(dat);
    if length(fn) == 3
        g = 1:3;
    else
        g = 1;
    end
    win = 180; %seconds (5 minutes)
    step = 10; %seconds
    avt = 50*step;
    nwin = 50*win;
    for ii = g
        if i == 1 || i == 4 %just HTA1 and VTA2_vp1
            bins = 1:5;
        else
            bins = 9:23;
        end
        x = dat.(fn{ii}).x;
        y = dat.(fn{ii}).y;
        nsamp = length(x);
        ind = [1 avt:avt:nsamp];
        for j = 1:length(ind)
            if abs(nsamp-ind(j)) < nwin  %skip the last few indexes approaching the end of the t-s
                continue
            else
                idx = ind(j):ind(j)+nwin-1;
            end
            time = dat.(fn{ii}).time(idx(1));
            xx = mean(x(idx,bins),2);
            yy = mean(y(idx,bins),2);
            rot = (pi*heading(i))/180;
            T = [cos(rot) -sin(rot);...
                sin(rot) cos(rot)];
            vels = [xx yy];
            V = vels*T';
            xx = V(:,1);yy = V(:,2);
            
            %calculate RMS velocities (Luhar et al. 2013)
            Ec = (1/length(yy))*sum(yy);Nc = (1/length(xx))*sum(xx);
            Ewrms = sqrt((1/length(yy))*sum((yy-Ec).^2));
            Nwrms = sqrt((1/length(xx))*sum((xx-Nc).^2));
            Uc = sqrt(Ec^2+Nc^2);Uwrms = sqrt(Ewrms^2+Nwrms^2);
            Uw = sqrt(2)*Uwrms;
            
            if i == 1 && ii == 3
                ff = fpeaks(1,:);
            elseif i == 2 && ii == 2
                ff = fpeaks(2,:);
            elseif i == 4
                ff = fpeaks(3,:);
            else
                ff = [NaN NaN NaN];
            end
            dd = d(i,ii);
            Re = (Uw*dd)/1.05E-6;
            St = (ff.*dd)./Uw;
            Lt = Uw./(2*pi*ff);
            %save variables
            data.(days{i}).(fn{ii}).time(j) = time;
            data.(days{i}).(fn{ii}).Uc(j) = Uc;
            data.(days{i}).(fn{ii}).Uw(j) = Uw;
            data.(days{i}).(fn{ii}).Re(j) = Re;
            data.(days{i}).(fn{ii}).d(j) = dd;
            data.(days{i}).(fn{ii}).St(j,:) = St;
            data.(days{i}).(fn{ii}).Lt(j,:) = Lt;
        end
    end
end
save('d:\Projects\Mekong_W2015\DataAnalysis\Paper2\HTAVTA_UwUcStLt','data','-v7.3')