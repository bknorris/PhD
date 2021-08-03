%Calculate bottom trace from VPs, then calculate bottom variance,
%normalized bottom variance (Staudt et al. 2017). Finally, try to correlate
%velocity measurements with bottom movement with a cross spectra (e.g.
%Puleo et al. 2014b).
%
%Paper #3

clear

fdir = 'd:\Projects\Mekong_W2015\Data\Vectrino\5March2015\';
files = dir([fdir '*March*']);files = {files.name};
vph = [0.125 0.125 0.08]; %vectrino height

for i = 1:length(files)
    disp(['Loading ' files{i}])
    load([fdir files{i}])
    
    %calculate bottom distance
    rb = VPRO.Data.Profiles_Range;
    h = vph(i)-rb;
    bd = VPRO.Data.BottomCheck_BottomDistance;
    bdt = VPRO.Data.BottomCheck_HostTimeMatlab;
    
    %run running median
    bdmid = my_running_median(bd,512);
    bdlvl = vph(i)-bdmid;
    bdlvl = smooth(bdlvl,10,'sgolay',3);
    
    %calc bottom variance
    n = 10*60*10;
    avt = n:n:length(bdlvl);
    win = [1 avt];
    bvar = zeros(length(win),1);
    for j = 1:length(win)-1
        bvar(j) = var(bdlvl(win(j):win(j+1)));
    end

    
        
        