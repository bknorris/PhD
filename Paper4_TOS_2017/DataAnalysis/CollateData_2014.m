clear
%Create a data file for the 2014 data (Aaron's, ours)
UWaqd = {'AQD_SW_2014.mat'};
WKaqd = {'HR3.mat';'HR4.mat'};
UWrbr = {'RBR_SW_2014.mat'};
UWhobo = 'SW_2014.mat';
data = struct();

dir1 = 'd:\Projects\Mekong_W2015\DataAnalysis\TOS\';
dir2 = 'd:\Projects\Mekong_F2014\Data\Aquadopp\F2F\';
%Load UW Aquadopps, filter above-surface vels,
%compute depth averages, compute SSC
%AQD Dep info:
lat = 9.494080; %I made these up based on Aaron's site description (150m along transect)
offset = 0.3150;
nn = {'sw'};
for i = 1:length(UWaqd)
    disp(['Loading ' dir1 UWaqd{i}])
    load([dir1 UWaqd{i}])
    
    %run surface tracking:
    
    %calc depth
    p = aqd.pres-offset(i);
    x = lat(i);
    g = zeros(length(p),1);h = zeros(length(p),1);
    for j = 1:length(p)
        g(j,:) = 9.780318*(1.0+(5.2788E-3+2.36E-5*x)*x)+1.092E-6*p(j,:);
        h(j,:) = ((((-1.82E-15*p(j,:)+2.279E-10)*p(j,:)-2.2512E-5)*p(j,:)+9.72659)*p(j,:))/g(j,:);
    end
    
    rb = repmat(aqd.z,length(h),1);
    for j = 1:length(h)
        idx = rb(j,:) >= h(j);
        aqd.v1(j,idx) = NaN;
        aqd.v2(j,idx) = NaN;
    end
    u = aqd.v1;v = aqd.v2;
    %calc velocity magnitude
    vel = sqrt(u.^2+v.^2);
    vel = nanmean(vel,2);
    time = aqd.time;
    
    %calc SSC (from Aaron)
    SSC = 0.0688.*aqd.ext2+50;
    SSC(SSC<450) = NaN;
    
    name = 'aqd';
    name2 = ['UW' nn{i}];
    %save data to structure
    data.(name).(name2).time = time;
    data.(name).(name2).vel = vel;
    data.(name).(name2).depth = h;
    data.(name).(name2).SSC = SSC;
    clear aqd
end

%Load WK Aquadopps, compute depth averages, compute depth, compute Hs
lat = [9.491309 9.491882];
zp = [0.092 0.082];zuv = zp;
nn = {'md';'fr'};
sr = 8;
for i = 1:length(WKaqd)
    name = 'aqd';
    name2 = ['WK' nn{i}];
    disp(['Loading ' dir2 WKaqd{i}])
    load([dir2 WKaqd{i}])
    time = aqdp.datenum;
    p = aqdp.pressure;
    u = aqdp.u;
    v = aqdp.v;

    %now average into 10-min chunks like Aaron's data
    avt = sr*300; %2 minute step
    nwin = sr*600; %10 minute window
    swin = sr*60; %60 second averaging window (to smooth)
    nsamp = length(p);
    ind = [1 avt:avt:nsamp];
    x = lat(i);
    lf = 1.2;
    hf = 0.05;
    for ii = 1:length(ind)
        if abs(nsamp-ind(ii)) < nwin  %skip the last few indexes approaching the end of the t-s
            continue
        else
            idx = ind(ii):ind(ii)+nwin-1;
        end
        %spectra of ADVs
        uu = u(idx,:);uu = cmgbridge(uu,100,1000,10000);
        vv = v(idx,:);vv = cmgbridge(vv,100,1000,10000);
        vels = sqrt(uu.^2+vv.^2);
        vels = nanmean(nanmean(vels,2));
        pp = p(idx)+zp(i);
        pp = cmgbridge(pp,100,1000,10000);
        tt = time(idx(length(idx)/2));
        
        g = zeros(length(pp),1);h = zeros(length(pp),1);
        for j = 1:length(pp)
            g(j,:) = 9.780318*(1.0+(5.2788E-3+2.36E-5*x)*x)+1.092E-6*pp(j,:);
            h(j,:) = ((((-1.82E-15*pp(j,:)+2.279E-10)*pp(j,:)-2.2512E-5)*pp(j,:)+9.72659)*pp(j,:))/g(j,:);
        end
        h = nanmean(h);
        g = nanmean(g);
        pp = detrend(pp);
        %sig wave height should be calculated with non-directional
        %spectra, i.e. from the pressure t-s
        [Cpp,F] = pwelch(pp,hanning(swin),swin*0.5,swin,sr);
        
        %wave parameters
        lfc = find(F >= hf,1,'first');hfc = find(F <= lf,1,'last');
        df = F(3)-F(2);
        omega = 2*pi.*F;
        k = qkhf(omega,h)./h;
        kh = k*h;
        kz = k*zuv(i);
        attn = cosh(kz)./cosh(kh);
        attn(attn<0.2) = 0.2;
        Spp = Cpp./(attn.^2);           %surface elevation spectrum
        m0 = sum(Spp(lfc:hfc)*df);
        Hs = 4*sqrt(m0);                %sig. wave height
        
        %save data to structure
        data.(name).(name2).time(ii) = tt;
        data.(name).(name2).vel(ii) = vels;
        data.(name).(name2).depth(ii) = h;
        data.(name).(name2).Hs(ii) = Hs;
        clear aqdp
    end
end

%Load UW RBRs, save depth, Hs, and calc SSC
nn = {'sw'};
for i = 1:length(UWrbr)
    name = 'rbr';
    name2 = ['UW' nn{i}];
    disp(['Loading ' dir1 UWrbr{i}])
    load([dir1 UWrbr{i}])
    
    %calc SSC
    SSC = 1.6252.*RBR_Data.turb-68.382;
    
    %save data to structure
    data.(name).(name2).time = RBR_Data.time;
    data.(name).(name2).SSC = SSC;
    data.(name).(name2).Hs = RBR_Data.Hsig;
    data.(name).(name2).depth = RBR_Data.depth;
    if i == 2
        data.(name).(name2).salt = RBR_Data.sal;
    end
end

save('d:\Projects\Mekong_W2015\DataAnalysis\TOS\Data_2014','data','-v7.3')