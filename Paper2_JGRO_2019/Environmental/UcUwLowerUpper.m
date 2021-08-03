%Calculate Uc and Uw (e.g., Luhar et al. 2013). For the HTA experiments, Uc
%and Uw come from the upper ADV (5109) for the 'above canopy' measurements
%and from VP1 for the 'within canopy' measurements. For the VTA, these come
%from the upper and lower VP, respectively. This is Version 1 of this
%script; to be used in tandem with the plot routine 'InAboveCanopyWhisks'

clear
%HTA instruments
inst{1} = 'D:\Projects\Mekong_W2015\Data\Vector\FSS\V5109_070315.mat';
inst{2} = 'D:\Projects\Mekong_W2015\Data\Vector\FSS\V5109_080315.mat';
inst{3} = 'D:\Projects\Mekong_W2015\Data\Vector\FSS\V5109_100315.mat';
inst{4} = 'D:\Projects\Mekong_W2015\DataAnalysis\Paper3\VPs\7March2015_Vels.mat';
inst{5} = 'D:\Projects\Mekong_W2015\DataAnalysis\Paper3\VPs\8March2015_Vels.mat';
inst{6} = 'D:\Projects\Mekong_W2015\DataAnalysis\Paper3\VPs\10March2015_Vels.mat';
%VTA instruments
inst{7} = 'D:\Projects\Mekong_W2015\Data\RBR\Duet\DPS2\Duet_140315.mat';
inst{8} = 'D:\Projects\Mekong_W2015\DataAnalysis\Paper3\VPs\14March2015a_Vels.mat';
exp = {'day1';'day2';'day3';'day4'};
%Initialize variables
for i = 1:3
    data.sw.(exp{i}) = struct();
end
data.ne.(exp{4}) = struct();
nwin = 20*60;
step = 10*60;
%HTA experiments - ADV first
for i = 1:3
    disp(['Loading: ' inst{i}])
    load(inst{i})
    time = ADV.datetime;
    u = ADV.U;
    v = ADV.V;
    p = ADV.Pres;
    %Downsample time series to 2 Hz
    u = navg(u,16);
    v = navg(v,16);
    p = navg(p,16);
    time = downsample(time,16);
    if length(time) > length(u) %time and velocity must be the same length.
        time = time(1:length(u));
    end
    zp = ADV.Metadata.pressure_sensor_height/1000; %in m
    zuv = ADV.Metadata.velocity_height/1000;
    lat = ADV.Metadata.inst_lat;
    %Run the full experiment (longest t-s); each file has two tides
    nsamp = length(time);
    ind = [1 step:step:nsamp];
    for ii = 1:length(ind)
        if abs(nsamp-ind(ii)) < nwin
            continue
        else
            idx = ind(ii):ind(ii)+nwin-1;
            upu = cmgbridge(u(idx),100,1000,1000);
            upv = cmgbridge(v(idx),100,1000,1000);
            upp = cmgbridge(p(idx),100,1000,1000);
            data.sw.(exp{i}).utime(ii) = time(ind(ii));
        end
        if nnz(isnan(upu))/length(upu) > 0.1 || nnz(isnan(upv))/length(upv) > 0.1 || nnz(isnan(upp))/length(upp) > 0.1
            data.sw.(exp{i}).uwspd(ii) = NaN;
            data.sw.(exp{i}).uwdir(ii) = NaN;
            data.sw.(exp{i}).Hs(ii) = NaN;
            data.sw.(exp{i}).Tr(ii) = NaN;
            data.sw.(exp{i}).h(ii) = NaN;
            data.sw.(exp{i}).ucspd(ii) = NaN;
            data.sw.(exp{i}).ucdir(ii) = NaN;
        else
            upu(isnan(upu)) = 0;
            upv(isnan(upv)) = 0;
            upp(isnan(upp)) = 0;
            [h,g] = estdepth(upp,lat);
            ws = puv(upp,upu,upv,mean(h),zp,zuv,2,round(nwin/8),1024,nwin,0.05,1.9);
            data.sw.(exp{i}).uwspd(ii) = ws.ubr;
            data.sw.(exp{i}).uwdir(ii) = wrapTo360(360-ws.azr);
            data.sw.(exp{i}).Hs(ii) = ws.Hrmsp*sqrt(2);
            data.sw.(exp{i}).Tp(ii) = ws.Tr;
            data.sw.(exp{i}).h(ii) = mean(h);
            s = uvwstats(upu,upv,[],0);
            data.sw.(exp{i}).ucspd(ii) = s.S;
            data.sw.(exp{i}).ucdir(ii) = s.az0;
        end
    end
    clear ADV
end
%HTA Experiments - VPs
loadings = 4:6;
for i = 1:3
    disp(['Loading: ' inst{loadings(i)}])
    m = matfile(inst{loadings(i)});
    fn = whos(m);fn = {fn.name};
    dat = m.(fn{2});
    time = dat.time;
    u = nanmean(dat.y,2); %VP x-beam points north
    v = nanmean(dat.x,2);
    %Downsample time series to 2 Hz
    u = navg(u,25);
    v = navg(v,25);
    time = downsample(time,25);
    if length(time) > length(u) %time and velocity must be the same length.
        time = time(1:length(u));
    end
    %Run the full experiment (longest t-s); measurements at beginning of
    %tide are likely a bit unreliable...
    nsamp = length(time);
    ind = [1 step:step:nsamp];
    for ii = 1:length(ind)
        if abs(nsamp-ind(ii)) < nwin
            continue
        else
            idx = ind(ii):ind(ii)+nwin-1;
            lou = cmgbridge(u(idx),100,1000,1000);
            lov = cmgbridge(v(idx),100,1000,1000);
            data.sw.(exp{i}).ltime(ii) = time(ind(ii));
        end
        if nnz(isnan(lou))/length(lou) > 0.1 || nnz(isnan(lov))/length(lov) > 0.1
            data.sw.(exp{i}).lwspd(ii) = NaN;
            data.sw.(exp{i}).lwdir(ii) = NaN;
            data.sw.(exp{i}).lcspd(ii) = NaN;
            data.sw.(exp{i}).lcdir(ii) = NaN;
        else
            lou(isnan(lou)) = 0;
            lov(isnan(lov)) = 0;
            s = uvwstats(lou,lov,[],0);
            data.sw.(exp{i}).lcspd(ii) = s.S;
            data.sw.(exp{i}).lcdir(ii) = s.az0;
        end
    end
    clear dat
end
%VTA experiments - VP3 and Duet
disp(['Loading: ' inst{7}])
load(inst{7})
disp(['Loading: ' inst{8}])
m = matfile(inst{8});
fn = whos(m);fn = {fn.name};
dat = m.(fn{3});
vpt = dat.time;
rbt = RBR.Datetime;
u = nanmean(dat.y,2); %VP x-beam points 20 deg east of N
v = nanmean(dat.x,2);
th = -20.*pi/180;
R = [cos(th) -sin(th); sin(th) cos(th)];
rxy = zeros(length(u),2);
for j = 1:length(u)
    rxy(j,:) = [v(j) u(j)]*R;
end
v = rxy(:,1);u = rxy(:,2);
p = RBR.SeaPres;
%Downsample time series to 2 Hz
u = navg(u,25);
v = navg(v,25);
vpt = downsample(vpt,25);
%Downsample time series to 2 Hz
p = navg(p,8);
rbt = downsample(rbt,8);
id = find(rbt>=vpt(1)&rbt<=vpt(end));
p = p(id);rbt = rbt(id);
time = (vpt(1):datenum(0,0,0,0,0,0.5):vpt(end))';
p = interp1(rbt,p,time);
if length(u) && length(v) ~= length(time) %#ok<ISMT>
    u = u(1:length(time));
    v = v(1:length(time));
end
zp = 0.52; %height of the duet
zuv = 0.806-0.04-0.015; %middle of VP sample volume
lat = 9.565495; %latitude of VTA
%Run the full experiment (longest t-s);
nsamp = length(time);
ind = [1 step:step:nsamp];
for ii = 1:length(ind)
    if abs(nsamp-ind(ii)) < nwin
        continue
    else
        idx = ind(ii):ind(ii)+nwin-1;
        upu = cmgbridge(u(idx),100,1000,1000);
        upv = cmgbridge(v(idx),100,1000,1000);
        upp = cmgbridge(p(idx),100,1000,1000);
        data.ne.(exp{4}).utime(ii) = time(ind(ii));
    end
    if nnz(isnan(upu))/length(upu) > 0.1 || nnz(isnan(upv))/length(upv) > 0.1 || nnz(isnan(upp))/length(upp) > 0.1
        data.ne.(exp{4}).uwspd(ii) = NaN;
        data.ne.(exp{4}).uwdir(ii) = NaN;
        data.ne.(exp{4}).Hs(ii) = NaN;
        data.ne.(exp{4}).Tr(ii) = NaN;
        data.ne.(exp{4}).h(ii) = NaN;
        data.ne.(exp{4}).ucspd(ii) = NaN;
        data.ne.(exp{4}).ucdir(ii) = NaN;
    else
        upu(isnan(upu)) = 0;
        upv(isnan(upv)) = 0;
        upp(isnan(upp)) = 0;
        [h,g] = estdepth(upp,lat);
        ws = puv(upp,upu,upv,mean(h),zp,zuv,2,round(nwin/8),1024,nwin,0.05,1.9);
        data.ne.(exp{4}).uwspd(ii) = ws.ubr;
        data.ne.(exp{4}).uwdir(ii) = wrapTo360(360-ws.azr);
        data.ne.(exp{4}).Hs(ii) = ws.Hrmsp*sqrt(2);
        data.ne.(exp{4}).Tp(ii) = ws.Tr;
        data.ne.(exp{4}).h(ii) = mean(h);
        s = uvwstats(upu,upv,[],0);
        data.ne.(exp{4}).ucspd(ii) = s.S;
        data.ne.(exp{4}).ucdir(ii) = s.az0;
    end
end
clear dat RBR
disp(['Loading: ' inst{8}])
dat = m.(fn{1});
time = dat.time;
u = nanmean(dat.y,2); %VP x-beam points 20 deg east of N
v = nanmean(dat.x,2);
rxy = zeros(length(u),2);
for j = 1:length(u)
    rxy(j,:) = [v(j) u(j)]*R;
end
v = rxy(:,1);u = rxy(:,2);
%Downsample time series to 2 Hz
u = navg(u,25);
v = navg(v,25);
time = downsample(time,25);
if length(time) > length(u) %time and velocity must be the same length.
    time = time(1:length(u));
end
%Run the full experiment (longest t-s); measurements at beginning of
%tide are likely a bit unreliable...
nsamp = length(time);
ind = [1 step:step:nsamp];
for ii = 1:length(ind)
    if abs(nsamp-ind(ii)) < nwin
        continue
    else
        idx = ind(ii):ind(ii)+nwin-1;
        lou = cmgbridge(u(idx),100,1000,1000);
        lov = cmgbridge(v(idx),100,1000,1000);
        data.ne.(exp{4}).ltime(ii) = time(ind(ii));
    end
    if nnz(isnan(lou))/length(lou) > 0.1 || nnz(isnan(lov))/length(lov) > 0.1
        data.ne.(exp{4}).lwspd(ii) = NaN;
        data.ne.(exp{4}).lwdir(ii) = NaN;
        data.ne.(exp{4}).lcspd(ii) = NaN;
        data.ne.(exp{4}).lcdir(ii) = NaN;
    else
        lou(isnan(lou)) = 0;
        lov(isnan(lov)) = 0;
        s = uvwstats(lou,lov,[],0);
        data.ne.(exp{4}).lcspd(ii) = s.S;
        data.ne.(exp{4}).lcdir(ii) = s.az0;
    end
end
clear dat
disp('Finished!')
save('d:\Projects\Mekong_W2015\DataAnalysis\Paper2\WaveStats_Full\UcUwUpLow_HTA_VTA','-struct','data','-v7.3')

