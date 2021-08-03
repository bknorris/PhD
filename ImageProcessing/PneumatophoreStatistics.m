%Pneumatophore statistics, calculates n, a, and phi, plus basal diameters
%(mean, max), mean and max canopy height

%28/03/2017
clear,close all
disp('Running vegetation statistics calculator')
disp('Select a working directory that contains vegetation statistics')
data = struct();
for fd = 1:2
    if fd == 1
        folder = 'd:\Mekong_F2014\DataAnalysis\Paper1\VegStats';
    elseif fd == 2
        folder = 'd:\Mekong_W2015\DataAnalysis\Paper1\VegStats';
    end
    
    % folder = uigetdir;
    fname = dir([folder '\' '*stats_fixed.mat']);fname = {fname.name};
    vegdat = struct();
    
    %Load the VP positions data
    vpdir = 'd:\Mekong_W2015\DataAnalysis\DataReports\';
    load([vpdir 'VP_positions.mat'])
    if strfind(folder,'2014') > 1
        VPpos.Qname = VPpos.Qname(1:21);
        VPpos.VP_X = VPpos.VP_X(1:21);
        VPpos.VP_Y = VPpos.VP_Y(1:21);
        VPpos.VP_Z = VPpos.VP_Z(1:21);
        year = 1;
    elseif strfind(folder,'2015') > 1
        VPpos.Qname = VPpos.Qname(22:end);load
        VPpos.VP_X = VPpos.VP_X(22:end);
        VPpos.VP_Y = VPpos.VP_Y(22:end);
        VPpos.VP_Z = VPpos.VP_Z(22:end);
        year = 2;
    end
    %%%Rotate quadrat so the top of fig points N
    if year == 1
        rotx = [0 0 180 0 0 0 0 0 0 0 0 0 0 0 180 0 0];
        roty = [0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 180 0];
        rotz = [90 90 0 -90 0 0 0 0 0 0 0 0 0 0 90 90 0];
    elseif year == 2
        rotx = [0 180 180 0 0 180 0 180];
        roty = [0 0 0 0 180 0 0 180];
        rotz = [90 0 0 270 90 0 -90 0];
    end
    for k = 1:length(fname)
        %%%Load the Data
        disp(['Loading ' folder '\' fname{k}])
        load([folder '\' fname{k}])
        qname = regexprep(fname{k},'stats_fixed.mat','');
        DATA = rmfield(DATA,'info'); %remove info field for processing
        
        %%%Load the correct VP data and scatterplot
        vpid = zeros(length(VPpos.Qname),1);
        for b = 1:length(VPpos.Qname)
            vpid(b) = strcmp(VPpos.Qname{b},qname);
        end
        vpx = (VPpos.VP_X.*vpid)';vpx(vpx == 0) = [];
        vpy = (VPpos.VP_Y.*vpid)';vpy(vpy == 0) = [];
        vpz = (VPpos.VP_Z.*vpid)';vpz(vpz == 0) = [];
        
        %%%First, plot the entire quadrat as an overview
        f1 = figure(1);
        set(f1,'PaperOrientation','portrait',...
            'position',[200 300   800   700]);
        
        %Create rotation matrix, rotate about z-axis if specified, flip
        %axes if specified
        thx = rotx(k)*pi/180;thy = roty(k)*pi/180;thz = rotz(k)*pi/180;
        Tx = [cos(2*thx) sin(2*thx); sin(2*thx) -cos(2*thx)];
        Ty = [-cos(2*thy) sin(2*thy); sin(2*thy) cos(2*thy)];
        Tz = [cos(thz) -sin(thz); sin(thz) cos(thz)];
        
        %Adjust VP positions based on rotations and reflections
        o = repmat([0.5;0.5],1,length(vpx)); %origin
        vpxy = [vpx;vpy];
        if rotz(k) ~= 0
            vpxy = Tz*(vpxy-o)+o;
        end
        if rotx(k) ~= 0
            vpxy = Tx*(vpxy-o)+o;
        end
        if roty(k) ~= 0
            vpxy = Ty*(vpxy-o)+o;
        end
        
        %Plot VP Positions
        pp = plot(vpxy(1,:),vpxy(2,:),'^','MarkerSize',12);
        set(pp,'MarkerEdgeColor','k',...
            'MarkerFaceColor',[1 0.9 0])
        hold on
        l = length(DATA.layer1.Rcenter);
        x = DATA.layer1.Rcenter(:,1)';
        y = DATA.layer1.Rcenter(:,2)';
        
        %Rotate & Reflect coordinates
        o = repmat([0.5;0.5],1,length(x));
        xy = [x;y];
        if rotz(k) ~= 0
            xy = Tz*(xy-o)+o;
        end
        if rotx(k) ~= 0
            xy = Tx*(xy-o)+o;
        end
        if roty(k) ~= 0
            xy = Ty*(xy-o)+o;
        end
        
        %Plot basal diameters
        for i = 1:l
            xc = xy(1,i);
            yc = xy(2,i);
            r = mean(DATA.layer1.Rradii(i,:));
            ang=0:0.01:2*pi;
            xp=r*cos(ang);
            yp=r*sin(ang);
            p = plot(xc+xp,yc+yp);set(p,'linewidth',1.5)
            hold on
        end
        axis equal
        grid on
        set(gca,'Xlim',[0 1],'Ylim',[0 1],'grid',':')
        xlabel('X'),ylabel('Y')
        title(qname)
        hold off
        
        %calculate pneumatophore statistics
        fn =fieldnames(DATA);
        v = length(fn);
        [vv,~] = size(DATA.layer1.Rcenter);
        Acir = zeros(v,vv);
        acir = zeros(v,1);
        D = zeros(v,vv);
        dmean = zeros(v,1);
        dmin = zeros(v,1);
        dmax = zeros(v,1);
        dstd = zeros(v,1);
        N = zeros(v,1);
        dist = zeros(v,vv);
        deltaS = zeros(v,1);
        for i = 1:v
            [n,~] = size(DATA.(fn{i}).Rradii);
            x = DATA.(fn{i}).Rcenter(:,1);
            y = DATA.(fn{i}).Rcenter(:,2);
            for j = 1:n
                r = mean(DATA.(fn{i}).Rradii(j,:));
                Acir(i,j) = pi*r^2;
                D(i,j) = 2*r;
                x1 = x(j);y1 = y(j);
                if length(x) > 1
                    xx = setxor(x,x1);yy = setxor(y,y1);
                    dd = sqrt((xx-x1).^2+(yy-y1).^2);
                    dist(i,j) = min(dd);
                else
                    dist(i,j) = 0;
                end
            end
            N(i,:) = n;
            zid = D == 0;
            D(zid) = NaN; %remove zeros, they compromise mean(D)
            zid = dist == 0;
            dist(zid) = NaN;
            deltaS(i,:) = nanmean(dist(i,:));
            dmean(i,:) = nanmean(D(i,:));
            dmin(i,:) = nanmin(D(i,:));
            dmax(i,:) = nanmax(D(i,:));
            dstd(i,:) = nanstd(D(i,:));
            acir(i,:) = sum(Acir(i,:));
        end
        zid = isnan(deltaS);deltaS(zid) = 0;
        zid = isnan(dmean);dmean(zid) = 0; %replace zeros
        zid = isnan(dmin);dmin(zid) = 0;
        zid = isnan(dmax);dmax(zid) = 0;
        zid = isnan(dstd);dstd(zid) = 0;
        h = 0.005;
        area = 1;                        %1 m^2
        a = (N/area).*dmean;		%24/09/2016 definition for N is "stems/area"
        phi = acir/area;
        Hmax = h*v; %doesn't include offset from bed of first height slice
        %calculate h_mean
        diff = zeros(v,1);
        H = zeros(v,vv);
        for i = 1:v-1
            df = N(i)-(N(i+1));
            midpts = (DATA.(fn{i}).Layermaxh+DATA.(fn{i}).Layerminh)/2;
            if df < 0
                H(i) = NaN;
                diff(i) = NaN;
            else
                if df == 1
                    H(i,1) = midpts;
                    diff(i) = df;
                else
                    H(i,1:df) = midpts;
                    diff(i) = df;
                end
            end
        end
        H(end,1) = Hmax-h; %include max canopy height 
        Hz = H(H~=0);Hn = Hz(~isnan(Hz));
        H = sort(Hn);
        z = h:h:Hmax;
        
        %save variables
        name = regexprep(fname{k},'stats_fixed.mat','');
        yr = {'y14','y15'};
        
        %for Julia only
        veg.n = N(1);
        veg.d = D(1,:)';
        veg.H = H;
        file = [name '_' yr{year}];
        save(['d:\Mekong_W2015\DataAnalysis\JuliasPaper\VegData\' file],'veg','-v7.3')
        
        data.(yr{year}).(name).z = z;
        data.(yr{year}).(name).Hmean = mean(H);
        data.(yr{year}).(name).sigmaH = std(H);
        data.(yr{year}).(name).Hmax = Hmax;
        data.(yr{year}).(name).N = N;
        data.(yr{year}).(name).dmean = dmean;
        data.(yr{year}).(name).dmin = dmin;
        data.(yr{year}).(name).dmax = dmax;
        data.(yr{year}).(name).dstd = dstd;
        data.(yr{year}).(name).a = a;
        data.(yr{year}).(name).phi = phi;
        data.(yr{year}).(name).deltaS = deltaS;
        pause(1)
    end
end
%plot routine, plot some example figures:
keep data yr
fpath = 'd:\Mekong_W2015\Figures\Quadrats\Statistics\';
close all
%d statistics
for i = 1:2
    fn = fieldnames(data.(yr{i}));
    n = length(fn);
    rng default
    r = randi([5 14],1,n); %generate numbers for plotting herrorbars
    c = brewermap(n,'paired');
    b = zeros(n,1);
    figure
    for ii = 1:n
        z = data.(yr{i}).(fn{ii}).z;
        dmin = data.(yr{i}).(fn{ii}).dmin;
        dmean = data.(yr{i}).(fn{ii}).dmean;
        dmax = data.(yr{i}).(fn{ii}).dmax;
        nn = length(dmean);
        stp = 1:r(ii):nn;
        b(ii) = plot(dmean,z,'linewidth',1.5,'color',c(ii,:));hold on
        p = herrorbar(dmean(stp),z(stp),dmin(stp),dmax(stp),'.');
        set(p,'color',c(ii,:))
    end
    leg = legend(b,fn);
    set(gca,'ylim',[0 1])
    set(gcf,'color','w')
    year = {'2014';'2015'};
    title(['{\itd}, Diameter stats, ' year{i} ' Quadrats'])
    ylabel('z')
    xlabel('d (m)')
end
%n statistics
for i = 1:2
    fn = fieldnames(data.(yr{i}));
    n = length(fn);
    c = brewermap(n,'paired');
    b = zeros(n,1);
    figure
    for ii = 1:n
        z = data.(yr{i}).(fn{ii}).z;
        N = data.(yr{i}).(fn{ii}).N;
        b(ii) = plot(N,z,'linewidth',1.5,'color',c(ii,:));hold on
    end
    leg = legend(b,fn);
    set(gca,'ylim',[0 1])
    set(gcf,'color','w')
    year = {'2014';'2015'};
    title(['{\itN}, No. of Stems, ' year{i} ' Quadrats'])
    ylabel('z')
    xlabel('N')
end
%a statistics
for i = 1:2
    fn = fieldnames(data.(yr{i}));
    n = length(fn);
    c = brewermap(n,'paired');
    b = zeros(n,1);
    figure
    for ii = 1:n
        z = data.(yr{i}).(fn{ii}).z;
        a = data.(yr{i}).(fn{ii}).a;
        b(ii) = plot(a,z,'linewidth',1.5,'color',c(ii,:));hold on
    end
    leg = legend(b,fn);
    set(gca,'ylim',[0 1])
    set(gcf,'color','w')
    year = {'2014';'2015'};
    title(['{\ita}, Frontal Area, ' year{i} ' Quadrats'])
    ylabel('z')
    xlabel('a (m^-^1)')
end
%phi statistics
for i = 1:2
    fn = fieldnames(data.(yr{i}));
    n = length(fn);
    c = brewermap(n,'paired');
    b = zeros(n,1);
    figure
    for ii = 1:n
        z = data.(yr{i}).(fn{ii}).z;
        a = data.(yr{i}).(fn{ii}).a;
        b(ii) = plot(a,z,'linewidth',1.5,'color',c(ii,:));hold on
    end
    leg = legend(b,fn);
    set(gca,'ylim',[0 1])
    set(gcf,'color','w')
    year = {'2014';'2015'};
    title(['{\it\phi}, Volume Fraction, ' year{i} ' Quadrats'])
    ylabel('z')
    xlabel('\phi')
end
tt = {'d';'d';'n';'n';'a';'a';'phi';'phi'};
year = {'2015';'2014';'2015';'2014';...
    '2015';'2014';'2015';'2014';};
handles = findobj('type','figure');
for i = 1:length(handles);
    id = handles(i);
    export_fig([fpath tt{i} year{id}],handles(id),'-jpeg','-nocrop')
end
save('d:\Mekong_W2015\DataAnalysis\DataReports\Vegetation\Vegdat_14_15_1m','data','-v7.3')