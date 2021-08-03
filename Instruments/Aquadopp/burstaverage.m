clear
dirc = 'c:\Users\bkn5\Projects\Mekong_F2014\Data\Aquadopp\F2F\';
files = dir([dirc 'HR4.mat']);
for i = 1:length(files)
    filelist = {files.name};varname = [filelist{i}(1:end-4)];
    load([dirc varname])
    %do the burst averaging
    spb = aqdp.metadata.spb;
    burstl = (length(aqdp.burst)/spb);
    burstc = ceil(burstl);
    remainder = burstl - burstc; %get the fraction of samples left off at the end
    rind1 = spb*burstc-(spb-1);
    rind2 = spb*burstl;
    aqdp.u(isnan(aqdp.u)) = 0;
    aqdp.v(isnan(aqdp.v)) = 0;
    aqdp.w(isnan(aqdp.w)) = 0;
    uav = [];
    vav = [];
    wav = [];
    for i = 1:burstc
        ix1 = spb*i-(spb-1);
        ix2 = spb*i; %this value must always be n+4096 > than ix1
        if i == max(burstc)
            ind = (rind1:rind2);
        else
            ind = (ix1:ix2);
        end
        uav(i,:) = (sum(aqdp.u(ind,:))./numel(aqdp.u(ind)));
        vav(i,:) = (sum(aqdp.v(ind,:))./numel(aqdp.v(ind)));
        wav(i,:) = (sum(aqdp.w(ind,:))./numel(aqdp.w(ind)));
        yearav(i,:) = (sum(aqdp.yearday(ind,:))./numel(aqdp.yearday(ind)));
    end
    aqdp.uav = uav;aqdp.vav = vav;aqdp.wav = wav;
    [a,b] = size(aqdp.uav);
    aqdp.nSamp = length(aqdp.burst);
    aqdp.nBursts = a;
    aqdp.yearav = yearav;
    save([dirc varname '_p'],'aqdp')
    
    imagesc(aqdp.uav')
    colorbar
    pause(5)
    close all
end
disp('Finished!')