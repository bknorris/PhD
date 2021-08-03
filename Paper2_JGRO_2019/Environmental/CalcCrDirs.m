%Calculate current directions for the HTA and VTA
%Lower currents come from the Vectrinos, (lets say VP3, because it was
%right underneath the ADV), and upper currents from the ADV for the HTA.
%Lower and upper velocities are both Vectrinos, VTA.
clear
inst{1} = 'D:\Projects\Mekong_W2015\Data\Vector\FSS\V5109_070315.mat';
inst{2} = 'D:\Projects\Mekong_W2015\DataAnalysis\Paper3\VPs\7March2015_Vels.mat';
inst{3} = 'D:\Projects\Mekong_W2015\Data\Vector\FSS\V5109_080315.mat';
inst{4} = 'D:\Projects\Mekong_W2015\DataAnalysis\Paper3\VPs\8March2015_Vels.mat';
inst{5} = 'D:\Projects\Mekong_W2015\Data\Vector\FSS\V5109_100315.mat';
inst{6} = 'D:\Projects\Mekong_W2015\DataAnalysis\Paper3\VPs\10March2015_Vels.mat';
inst{7} = 'D:\Projects\Mekong_W2015\DataAnalysis\Paper3\VPs\14March2015a_Vels.mat';
loading = [1 2;3 4;5 6];
exp = {'day1';'day2';'day3';'day4'};
%Initialize variables
for i = 1:3
    cur.sw.(exp{i}) = struct();
end
cur.ne.(exp{4}) = struct();

start(1) = datenum(2015,03,07,13,36,32);stop(1) = datenum(2015,03,07,17,05,02);
start(2) = datenum(2015,03,08,14,14,22);stop(2) = datenum(2015,03,08,19,01,21);
start(3) = datenum(2015,03,10,14,45,06);stop(3) = datenum(2015,03,10,16,36,36);
%HTA experiments
for i = 1:3
    %Load the ADV
    load(inst{loading(i,1)})
    %Load only the third Vectrino
    m = matfile(inst{loading(i,2)});
    fn = whos(m);fn = {fn.name};
    dat = m.(fn{3});
    %Need to get the ADV and Vectrino on the same timebase...
    idx = find(ADV.datetime>=start(i)&ADV.datetime<=stop(i));
    uu = ADV.U(idx);vu = ADV.V(idx);
    adt = ADV.datetime(idx);
    idx = find(dat.time>=start(i)&dat.time<=stop(i));
    ul = nanmean(dat.y(idx,:),2);vl = nanmean(dat.x(idx,:),2);
    vpt = dat.time(idx);
    time = start(i):datenum(0,0,0,0,0,1):stop(i);time = time(1:end-1);
    %Downsample to 1 Hz
    uu = navg(uu,32);vu = navg(vu,32);
    ul = navg(ul,50);vl = navg(vl,50);
    %Average over 10-minute blocks @ 5-minute step
    nwin = 60*10;
    step = 60*5;
    nsamp = length(time);
    ind = [1 step:step:nsamp];
    for ii = 1:length(ind)
        if abs(nsamp-ind(ii)) < nwin
            continue
        else
            idx = ind(ii):ind(ii)+nwin-1;
            upu = cmgbridge(uu(idx),100,1000,1000);
            upv = cmgbridge(vu(idx),100,1000,1000);
            lou = cmgbridge(ul(idx),100,1000,1000);
            lov = cmgbridge(vl(idx),100,1000,1000);
            cur.sw.(exp{i}).time(ii) = time(ind(ii));
        end
        if nnz(isnan(upu))/length(upu) > 0.1 || nnz(isnan(upv))/length(upv) > 0.1
            cur.sw.(exp{i}).uspd(ii) = NaN;
            cur.sw.(exp{i}).udir(ii) = NaN;
        else
            upu(isnan(upu)) = 0;
            upv(isnan(upv)) = 0;
            up = cmgpca(upu,upv,[],0);
            cur.sw.(exp{i}).uspd(ii) = up.mspd;
            cur.sw.(exp{i}).udir(ii) = up.mdir;
        end
        if nnz(isnan(lou))/length(lou) > 0.1 || nnz(isnan(lov))/length(lov) > 0.1
            cur.sw.(exp{i}).lspd(ii) = NaN;
            cur.sw.(exp{i}).ldir(ii) = NaN;
        else
            lou(isnan(lou)) = 0;
            lov(isnan(lov)) = 0;
            low = cmgpca(lou,lov,[],0);
            cur.sw.(exp{i}).lspd(ii) = low.mspd;
            cur.sw.(exp{i}).ldir(ii) = low.mdir;
        end
    end
    clear ADV dat
end
%VTA experiments
m = matfile(inst{7});
fn = whos(m);fn = {fn.name};
vp1 = m.(fn{1});vp3 = m.(fn{3});
uu = nanmean(vp3.y,2);vu = nanmean(vp3.x,2);
ul = nanmean(vp1.y,2);vl = nanmean(vp1.x,2);
%Downsample to 1 Hz
uu = navg(uu,50);vu = navg(vu,50);
ul = navg(ul,50);vl = navg(vl,50);
%Average over 10-minute blocks @ 5-minute step
nwin = 60*10;
step = 60*5;
time = vp1.time(1):datenum(0,0,0,0,0,1):vp1.time(end);
nsamp = length(time);
ind = [1 step:step:nsamp];
for ii = 1:length(ind)
    if abs(nsamp-ind(ii)) < nwin
        continue
    else
        idx = ind(ii):ind(ii)+nwin-1;
        upu = cmgbridge(uu(idx),100,1000,1000);
        upv = cmgbridge(vu(idx),100,1000,1000);
        lou = cmgbridge(ul(idx),100,1000,1000);
        lov = cmgbridge(vl(idx),100,1000,1000);
        cur.ne.(exp{4}).time(ii) = time(ind(ii));
    end
    if nnz(isnan(upu))/length(upu) > 0.1 || nnz(isnan(upv))/length(upv) > 0.1
        cur.ne.(exp{4}).uspd(ii) = NaN;
        cur.ne.(exp{4}).udir(ii) = NaN;
    else
        upu(isnan(upu)) = 0;
        upv(isnan(upv)) = 0;
        up = cmgpca(upu,upv,[],0);
        cur.ne.(exp{4}).uspd(ii) = up.mspd;
        cur.ne.(exp{4}).udir(ii) = up.mdir;
    end
    if nnz(isnan(lou))/length(lou) > 0.1 || nnz(isnan(lov))/length(lov) > 0.1
        cur.ne.(exp{4}).lspd(ii) = NaN;
        cur.ne.(exp{4}).ldir(ii) = NaN;
    else
        lou(isnan(lou)) = 0;
        lov(isnan(lov)) = 0;
        low = cmgpca(lou,lov,[],0);
        cur.ne.(exp{4}).lspd(ii) = low.mspd;
        cur.ne.(exp{4}).ldir(ii) = low.mdir;
    end
end
clear vp1 vp3
disp('Finished!')
save('d:\Projects\Mekong_W2015\DataAnalysis\Paper2\WaveStats_Full\CurStats_HTA_VTAv2','-struct','cur','-v7.3')


