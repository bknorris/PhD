%Collate data for the TOS data figure 1
%Requirements: Depth, Depth-averaged-velocity, Hs, SSC, and Salinity
%Hs needs to be calculated individually for every Aquadopp
%SSC needs to be calculated individually for every instrument
clear

%data directories


dir1 = 'd:\Projects\Mekong_W2015\DataAnalysis\TOS\';
if ~exist([dir1 'Data_Fig1.mat'],'file')
    UWaqd = {'AQD_SW_2015.mat';'AQD_NE_2015.mat'};
    WKaqd = {'AD5116_9March2015.mat';'AD5116_12March2015.mat';...
        'AD5116_15March2015.mat'};
    WKadv = {'V5108_030315.mat';'V5108_080315.mat';...
        'V5108_120315.mat';'V5109_130315.mat'};
    UWrbr = {'RBR_SW_2015.mat';'RBR_NE_2015.mat'};
    WKrbr = 'RBR_15244.mat';
    UWhobo = 'SW_2015.mat';
    load('d:\Projects\Mekong_W2015\Data\Weather\Mekong2015_metstation.mat')
    data = struct();
    
    %Load UW Aquadopps, filter above-surface vels,
    %compute depth averages, compute SSC
    %AQD Dep info:
    lat = [9.492706,9.565666]; %I made these up based on Aaron's site description (150m along transect)
    offset = [0.3150,0.362];
    nn = {'sw';'ne'};
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
        if i == 1
            SSC = 0.0688.*aqd.ext1+50;
        else
            SSC = 0.0688.*aqd.ext2+50;
        end
        
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
    lat = [9.491844,9.565448];
    zp = [0.075,0.105];zuv = zp;
    nn = {'sw';'ne'};
    sr = 8;
    for i = 1:2
        name = 'aqd';
        name2 = ['WK' nn{i}];
        if i == 1
            disp(['Loading ' dir1 WKaqd{i}])
            load([dir1 WKaqd{i}])
            
            %concatenate F2F2 and F2F3 aqdp files together
            t1 = datenum(2015,03,05,12,00,00);
            t2 = datenum(2015,03,12,12,00,00);
            stp = datenum(0,0,0,0,0,1/sr);
            time = t1:stp:t2;
            [m,n] = size(aqdp.u);
            q = length(time);
            u = NaN(q,n);
            v = NaN(q,n);
            p = NaN(q,1);
            
            u(1:m,:) = aqdp.u;
            v(1:m,:) = aqdp.v;
            p(1:m) = aqdp.pressure;
            clear aqdp
            
            disp(['Loading ' dir1 WKaqd{i+1}])
            load([dir1 WKaqd{i+1}])
            [m,n] = size(aqdp.u);
            m = m-1;
            u(end-m:end,1:n) = aqdp.u;
            v(end-m:end,1:n) = aqdp.v;
            p(end-m:end) = aqdp.pressure;
        else
            disp(['Loading ' dir1 WKaqd{i+1}])
            load([dir1 WKaqd{i+1}])
            time = aqdp.datenum;
            p = aqdp.pressure;
            u = aqdp.u;
            v = aqdp.v;
        end
        
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
    
    %Load WK ADVs, compute depth averages, compute depth, compute Hs
    lat = [9.491583,9.565383];
    zp = [0.196,0.12];
    zuv = [0.315,0.29];
    nn = {'sw';'ne'};
    sr = 32;
    for i = 1:2
        name = 'adv';
        name2 = ['WK' nn{i}];
        if i == 1
            disp(['Loading ' dir1 WKadv{i}])
            load([dir1 WKadv{i}])
            
            %concatenate all files together
            t1 = datenum(2015,03,03,21,30,00);
            t2 = datenum(2015,03,15,06,25,22);
            stp = datenum(0,0,0,0,0,1/sr);
            time = t1:stp:t2;
            [m,n] = size(ADV.V);
            q = length(time);
            u = NaN(q,n);
            v = NaN(q,n);
            p = NaN(q,1);
            
            u(1:m,:) = ADV.U;
            v(1:m,:) = ADV.V;
            p(1:m) = ADV.Pres;
            clear ADV
            
            disp(['Loading ' dir1 WKadv{i+1}])
            load([dir1 WKadv{i+1}])
            tmp = abs(time-ADV.datetime(1));
            [~,id] = min(tmp);
            
            [m,n] = size(ADV.U);
            m = m-1;
            u(id:id+m) = ADV.U;
            v(id:id+m) = ADV.V;
            p(id:id+m) = ADV.Pres;
            clear ADV
            
            disp(['Loading ' dir1 WKadv{i+2}])
            load([dir1 WKadv{i+2}])
            [m,n] = size(ADV.U);
            m = m-1;
            u(end-m:end,1:n) = ADV.U;
            v(end-m:end,1:n) = ADV.V;
            p(end-m:end) = ADV.Pres;
            clear ADV
        else
            disp(['Loading ' dir1 WKadv{i+2}])
            load([dir1 WKadv{i+2}])
            time = ADV.datetime;
            p = ADV.Pres;
            u = ADV.U;
            v = ADV.V;
        end
        
        %now average into 10-min chunks like Aaron's data
        avt = sr*300; %2 minute step
        nwin = sr*600; %10 minute window
        swin = sr*20; %20 second averaging window (to smooth)
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
    nn = {'sw';'ne'};
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
    
    
    %Load WK RBRs, save salinity
    name = 'rbr';
    name2 = 'WKsw';
    disp(['Loading ' dir1 WKrbr])
    load([dir1 WKrbr])
    
    %save data to structure
    data.(name).(name2).time = RBR.Datetime;
    data.(name).(name2).salt = RBR.Salinity;
    
    %Load UW Hobo, save salinity
    name = 'hobo';
    name2 = 'UWsw';
    disp(['Loading ' dir1 UWhobo])
    load([dir1 UWhobo])
    
    %save data to structure
    data.(name).(name2).time = hobo.Datetime;
    data.(name).(name2).salt = hobo.Salinity;
    
    save('d:\Projects\Mekong_W2015\DataAnalysis\TOS\Data_Fig1','data','-v7.3')
else
    load([dir1 'Data_Fig1.mat'])
end

%%%Plot Routine%%%
f1 = figure(1);
set(f1,'PaperOrientation','portrait',...
    'position',[500 200 1000 600]);
c = [207 176 126;
    60 166 74;
    4 76 41]./255;
%Southwest:
%plot depth
sp(1) = subplot(521);
id1 = data.adv.WKsw.depth < 0.2;
data.adv.WKsw.depth(id1) = NaN;
plot(data.adv.WKsw.time,data.adv.WKsw.depth,...
    'linewidth',1.5,'color',(c(1,:))), hold on
id2 = data.aqd.WKsw.depth < 0.2;
data.aqd.WKsw.depth(id2) = NaN;
plot(data.aqd.WKsw.time,data.aqd.WKsw.depth,...
    'linewidth',1.5,'color',(c(2,:)))
id3 = data.aqd.UWsw.depth < 0.26;
data.aqd.UWsw.depth(id3) = NaN;
plot(data.aqd.UWsw.time,data.aqd.UWsw.depth,...
    'linewidth',1.5,'color',(c(3,:)))
text(0.02,0.8,'Depth (m)','Units','normalized')
%plot depth-av-vel
sp(3) = subplot(523);
data.adv.WKsw.vel(id1) = NaN;
plot(data.adv.WKsw.time,data.adv.WKsw.vel,...
    'linewidth',1.5,'color',(c(1,:))), hold on
data.aqd.WKsw.vel(id2) = NaN;
plot(data.aqd.WKsw.time,data.aqd.WKsw.vel,...
    'linewidth',1.5,'color',(c(2,:)))
data.aqd.UWsw.vel(id3) = NaN;
plot(data.aqd.UWsw.time,data.aqd.UWsw.vel,...
    'linewidth',1.5,'color',(c(3,:)))
text(0.02,0.8,'$\langle{\bar{u}}\rangle \quad{(m s^{-1})}$',...
    'interpreter','latex','Units','normalized')
%plot Hs
sp(5) = subplot(525);
id = data.adv.WKsw.Hs < 0.0025;
data.adv.WKsw.Hs(id) = NaN;
plot(data.adv.WKsw.time,data.adv.WKsw.Hs,...
    'linewidth',1.5,'color',(c(1,:))), hold on
plot(data.rbr.UWsw.time,data.rbr.UWsw.Hs,...
    'linewidth',1.5,'color',(c(2,:)))
text(0.02,0.8,'H_s (m)','Units','normalized')
%plot SSC
sp(7) = subplot(527);
plot(data.rbr.UWsw.time,data.rbr.UWsw.SSC,...
    'linewidth',1.5,'color',(c(2,:))), hold on
id = data.aqd.UWsw.SSC < 100;
data.aqd.UWsw.SSC(id) = NaN;
plot(data.aqd.UWsw.time,data.aqd.UWsw.SSC,...
    'linewidth',1.5,'color',(c(3,:)))
text(0.02,0.8,'SSC (mg l^-^1)','Units','normalized')
%plot Salinity
sp(9) = subplot(529);
id = data.rbr.WKsw.salt < 16;
data.rbr.WKsw.salt(id) = NaN;
plot(data.rbr.WKsw.time,data.rbr.WKsw.salt,...
    'linewidth',1.5,'color',(c(1,:))),hold on
id = data.hobo.UWsw.salt < 13;
data.hobo.UWsw.salt(id) = NaN;
plot(data.hobo.UWsw.time,data.hobo.UWsw.salt,...
    'linewidth',1.5,'color',(c(2,:)))
text(0.02,0.8,'Salinity (PSU)','Units','normalized')

%Northeast:
%plot depth
sp(2) = subplot(522);
id = data.adv.WKne.depth < 0.13;
data.adv.WKne.depth(id) = NaN;
pp(1) = plot(data.adv.WKne.time,data.adv.WKne.depth,...
    'linewidth',1.5,'color',(c(1,:))); hold on
id = data.rbr.UWne.depth < 0.122;
data.rbr.UWne.depth(id) = NaN;
pp(2) = plot(data.rbr.UWne.time,data.rbr.UWne.depth,...
    'linewidth',1.5,'color',(c(2,:)));
id2 = data.aqd.UWne.depth < 0.43;
data.aqd.UWne.depth(id2) = NaN;
pp(3) = plot(data.aqd.UWne.time,data.aqd.UWne.depth,...
    'linewidth',1.5,'color',(c(3,:)));
leg = legend(pp,{'Mudflat';'Fringe';'Forest'});
%plot depth-av-vel
sp(4) = subplot(524);
plot(data.adv.WKne.time,data.adv.WKne.vel,...
    'linewidth',1.5,'color',(c(1,:))), hold on
plot(data.aqd.WKne.time,data.aqd.WKne.vel,...
    'linewidth',1.5,'color',(c(2,:)))
plot(data.aqd.UWne.time,data.aqd.UWne.vel,...
    'linewidth',1.5,'color',(c(3,:)))
%plot Hs
sp(6) = subplot(526);
id = data.adv.WKne.Hs < 0.0045;
data.adv.WKne.Hs(id) = NaN;
plot(data.adv.WKne.time,data.adv.WKne.Hs,...
    'linewidth',1.5,'color',(c(1,:))), hold on
plot(data.rbr.UWne.time,data.rbr.UWne.Hs,...
    'linewidth',1.5,'color',(c(2,:)))
%plot SSC
sp(8) = subplot(528);
plot(data.rbr.UWne.time,data.rbr.UWne.SSC,...
    'linewidth',1.5,'color',(c(2,:))), hold on
id = find(data.aqd.UWne.SSC < 100);
data.aqd.UWne.SSC(id) = NaN;
plot(data.aqd.UWne.time,data.aqd.UWne.SSC,...
    'linewidth',1.5,'color',(c(3,:)))
%plot Salinity
sp(10) = subplot(5,2,10);
plot(data.rbr.UWne.time,data.rbr.UWne.salt,...
    'linewidth',1.5,'color',(c(2,:))),hold on

%global plot adjustments:

%axis lims, ticks
t1 = datenum(2015,03,03,21,30,00);
t2 = datenum(2015,03,15,06,25,22);
stp = datenum(0,0,1,0,0,0);
set(sp,'xtick',t1:stp:t2,'xlim',[t1 t2])
set([sp(2) sp(4) sp(6) sp(8) sp(10)],...
    'yticklabel',[])
set([sp(1) sp(2) sp(3) sp(4) sp(5) sp(6) sp(7) sp(8)],...
    'xticklabel',[])
set([sp(1) sp(2)],'ylim',[0 3])
set([sp(3) sp(4)],'ylim',[0 0.75],'ytick',0:0.25:0.5)
set([sp(5) sp(6)],'ylim',[0 1],'ytick',0:0.33:0.66)
set([sp(7) sp(8)],'ylim',[0 3000],'ytick',0:1000:2500)
set([sp(9) sp(10)],'ylim',[10 25],'ytick',10:5:25)

%positioning
set(sp(1),'position',[0.06 0.79 0.4 0.15])
set(sp(2),'position',[0.48 0.79 0.4 0.15])
set(sp(3),'position',[0.06 0.62 0.4 0.15])
set(sp(4),'position',[0.48 0.62 0.4 0.15])
set(sp(5),'position',[0.06 0.45 0.4 0.15])
set(sp(6),'position',[0.48 0.45 0.4 0.15])
set(sp(7),'position',[0.06 0.28 0.4 0.15])
set(sp(8),'position',[0.48 0.28 0.4 0.15])
set(sp(9),'position',[0.06 0.11 0.4 0.15])
set(sp(10),'position',[0.48 0.11 0.4 0.15])
set(leg,'position',[0.91 0.84 0.05 0.05])

%labeling
title(sp(1),'Southwest')
title(sp(2),'Northeast')
xlabel(sp(9),'Day in March, 2015')
xlabel(sp(10),'Day in March, 2015')
datetick(sp(9),'x','dd','keepticks','keeplimits')
datetick(sp(10),'x','dd','keepticks','keeplimits')


% finishing
% prettyfigures('text',13,'labels',14,'box',1,'tickdir','in')
% export_fig('d:\Projects\Documents\Writing\TOSpaper\Figures\DataTimeseries\DataTimeSeries_v3','-pdf','-nocrop')






