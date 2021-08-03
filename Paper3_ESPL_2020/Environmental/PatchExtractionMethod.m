%plot a figure showing the 'big D versus little d' estimation method,
%including a current+wave rose showing dominant current direction.
clear, close all
quads = {'Q3Bstats_fixed.mat';'Q2Astats_fixed.mat'};
aqfil = {'F2F3_11_AD5117.mat';'FSS_08_AD5116.mat'};
vegdir = 'd:\Mekong_W2015\DataAnalysis\DataReports\Vegetation\2015\';
aqdir = {'d:\Mekong_W2015\DataAnalysis\Paper3\WaveStats\11-03-15\';...
    'd:\Mekong_W2015\DataAnalysis\Paper3\WaveStats\08-03-15\'};
sdir = 'e:\GradSchool\DataAnalysis\Paper3\Figures\';
%Load the VP positions data
vpdir = 'd:\Mekong_W2015\DataAnalysis\DataReports\';
load([vpdir 'VP_positions.mat'])
VPpos.Qname = VPpos.Qname(22:end);
VPpos.VP_X = VPpos.VP_X(22:end);
VPpos.VP_Y = VPpos.VP_Y(22:end);
VPpos.VP_Z = VPpos.VP_Z(22:end);
year = 2;
wvp = 2;
%%%Rotate quadrat so the top of fig points N 
rotx = [180 180];
roty = [0 0];
rotz = [0 0];
bs = 0.2; %box size (20cm^2)
%%%Initialize the Figure
f1 = figure(1);
set(f1,'PaperOrientation','portrait',...
    'position',[400 200   800   550]);
set(gcf,'color','w','paperpositionmode','auto',...
    'renderer','painters')
vplots = [2 4];
aplots = [1 3];
sp = zeros(1,4);
for k = 1:length(quads)
    %% Load the Data
    disp(['Loading ' vegdir quads{k}])
    disp(['Loading ' aqdir{k} aqfil{k}])
    load([vegdir quads{k}])
    load([aqdir{k} aqfil{k}])
    qname = regexprep(quads{k},'stats_fixed.mat','');
    DATA = rmfield(DATA,'info'); %remove info field for processing
    
    %Plot AQDP wave rose
    sp(aplots(k)) = subplot(2,2,aplots(k));
    cmap = brewermap(100,'PuBu');
    Dir1 = 360-wave.Dir;
    h = WindRose(mean(Dir1,2),wave.ubr',...
        'anglenorth',0,...
        'angleeast',90,...
        'titlestring',[],...
        'legendtype',1,...
        'axes',sp(aplots(k)));

    %Load the correct VP data and scatterplot
    vpid = zeros(length(VPpos.Qname),1);
    for b = 1:length(VPpos.Qname)
        vpid(b) = strcmp(VPpos.Qname{b},qname);
    end
    vpx = (VPpos.VP_X.*vpid)';vpx(vpx == 0) = [];
    vpy = (VPpos.VP_Y.*vpid)';vpy(vpy == 0) = [];
    vpz = (VPpos.VP_Z.*vpid)';vpz(vpz == 0) = [];
    %Check to see if the selection is correct with the 'while' statement
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
    
    %Draw a rectangle around the desired points (for help: 'help rbbox')
    if length(vpxy) == 2
        vcx = vpxy(1);
        vcy = vpxy(2);
    else
        vcx = vpxy(1,wvp);
        vcy = vpxy(2,wvp);
    end
    inc = bs/2;
    x = [vcx+inc vcx-inc vcx-inc vcx+inc vcx+inc];
    y = [vcy+inc vcy+inc vcy-inc vcy-inc vcy+inc];
    x = unique(x);xlo = x(1);xhi = x(2);
    y = unique(y);ylo = y(1);yhi = y(2);
    
    %Crop DATA layer1 based on answer inputs, then replot
    xid = find(xy(1,:) >= xlo & xy(1,:) <= xhi);
    yid = find(xy(2,:) >= ylo & xy(2,:) <= yhi);
    id = intersect(xid,yid);
    nstem = length(id);
    
    sp(vplots(k)) = subplot(2,2,vplots(k));
    %Plot VP Positions
    pp = plot(vpxy(1,:),vpxy(2,:),'^','MarkerSize',12);
    set(pp,'MarkerEdgeColor','k',...
        'MarkerFaceColor',[1 0.9 0])
    hold on
    l = length(id);
    R = zeros(l,1);
    X = zeros(l,1);
    Y = zeros(l,1);
    for i = 1:l
        xc = xy(1,id(i));
        yc = xy(2,id(i));
        r = mean(DATA.layer1.Rradii(id(i),:));
        ang=0:0.01:2*pi;
        xp=r*cos(ang);
        yp=r*sin(ang);
        p = plot(xc+xp,yc+yp);set(p,'linewidth',1.5,'color','k')
        hold on
        R(i,:) = r;
        X(i,:) = xc;
        Y(i,:) = yc;
    end
    axis equal
    grid on
    set(gca,'Xlim',[xlo xhi],'Ylim',[ylo yhi])
    xlabel(sprintf('Across-shore\nDistance [m]'))
    ylabel(sprintf('Along-shore\nDistance [m]'))
end
prettyfigures('text',13,'labels',14,'box',1)
export_fig([sdir 'PatchExtractionMethod'],'-pdf','-painters','-nocrop')