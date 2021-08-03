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
heading = [20 20 20]; %headings for 2015 are 20 deg
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
    bdlvl = smooth(bdlvl,3);
    
    %rotate vels to xyz, then rotate to along and across shore
    [VPRO.Data,VPRO.Config] = VPro_coordinateTransform(VPRO.Data,VPRO.Config,'bx');
    x = VPRO.Data.Profiles_VelX;
    y = VPRO.Data.Profiles_VelY;
    rot = (pi*heading(i))/180;
    across = x.*(ones(size(x))*cos(rot)) + ...
        y.*(ones(size(y))*sin(rot));
    along = -y.*(ones(size(y))*sin(rot)) + ...
        x.*(ones(size(x))*cos(rot));
    
    %calculate RMS velocities (Luhar et al. 2013)
    disp('Calculating RMS velocities (Uc & Uw)')
    n = 5; %downsample 50Hz to 10Hz
    avt = n:n:length(along);
    win = [1 avt];lw = length(win);
    Uc = zeros(lw,1);
    Uw = zeros(lw,1);
    Umag = zeros(lw,1);
    for j = 1:lw-1
        idx = win(j):win(j+1);
        al = nanmean(along(idx),1);ac = nanmean(across(idx),1);
        Ec = (1/length(al))*sum(al);
        Nc = (1/length(ac))*sum(ac);
        Ewrms = sqrt((1/length(al))*sum((al-Ec).^2));
        Nwrms = sqrt((1/length(ac))*sum((ac-Nc).^2));
        Uc(j) = sqrt(Ec^2+Nc^2);Uwrms = sqrt(Ewrms^2+Nwrms^2);
        Uw(j) = sqrt(2)*Uwrms;
        Umag(j) = nanmean(sqrt(ac.^2+al.^2));
    end
    
    %calc bottom variance
    n = 10*60*10; %10-minute window
    avt = n:n:length(bdlvl);
    win = [1 avt];
    bvar = zeros(length(win),1);
    bnorm = zeros(length(win),1);
    for j = 1:length(win)-1
        bvar(j) = var(bdlvl(win(j):win(j+1)));
        bnorm(j) = bvar(j)./(nanmean(Umag(win(j):win(j+1)))*600);
    end
    
        
        