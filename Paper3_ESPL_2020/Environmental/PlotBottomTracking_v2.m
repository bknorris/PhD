%Plot bed movement during the HTA experiments (A re-do of the original
%'scour.m' plot
clear
figdir = 'd:\Projects\Mekong_W2015\Figures\Paper1\';
vph = [0.247 0.271 0.229];
datdir = 'D:\Projects\Mekong_W2015\Data\Vectrino\9March2015\';
fname = dir([datdir '*2015_*.mat']);
inst = cell(3,1);
for i = 1:3
    inst{i} = [datdir char({fname(i).name})];
end
c = jet(3);
sdir = 'd:\Projects\Mekong_W2015\DataAnalysis\Paper3\BottomTrack\09-03-15\';
if exist([sdir 'BDtrace_floodebb.mat'],'file') > 0
    disp('Loading BDtrace_floodebb.mat')
    load([sdir 'BDtrace_floodebb.mat'])
else
    for i = 1:3
        start = datenum(2015,03,09,02,10,00);stop = datenum(2015,03,09,06,50,00);
        vname = lower(inst{i}(62:64));
        
        disp(['Loading ' inst{i}])
        load(inst{i});
        bd = VPRO.Data.BottomCheck_BottomDistance;
        time = VPRO.Data.BottomCheck_HostTimeMatlab;
        id = find(time >= start & time <= stop);
        bd = bd(id);time = time(id);
        clear VPRO
        %remove zeros
        time(bd == 0) = [];
        bd(bd == 0) = [];
        bdmid = my_running_median(bd,512);
        %window to smooth for plotting
        win = 10*50;
        step = 10*2;
        idx = [1 step:step:length(bdmid)];
        bd2 = zeros(length(idx),1);
        time2 = zeros(length(idx),1);
        for ii = 1:length(idx)-1
            if abs(length(bdmid)-idx(ii)) < win
                continue
            else
                bwin = bdmid(idx(ii):idx(ii)+win);
                if i == 1
                    bd2(ii,:) = max(bwin);
                elseif i == 2
                    bd2(ii,:) = mean(bwin);
                elseif i == 3
                    bd2(ii,:) = min(bwin);
                end
                time2(ii,:) = time(idx(ii));
            end
        end
        bd2(bd2 == 0) = []; %remove trailing zeros
        time2(time2 == 0) = [];
       
        bdavs = runningmean(bd2,100);
        %CHECK
%         figure
%         plot(time,bd),hold on
%         plot(time,bdmid,'r')
%         plot(time2,bd2,'g')
%         plot(time2,bdavs,'m','linewidth',3)
%         
%         figure(4)
%         plot(time2,vph(i)-bdavs,'color',c(i,:)), hold on
        
        data.(vname).time = time2;
        data.(vname).bd = bdavs;
        data.(vname).ble = vph(i)-bdavs;
        clear VPRO
    end
save([sdir 'BDtrace_floodebb'],'data','-v7.3')
end

%%%Plot Routine%%%
f1 = figure(1);
set(f1,'PaperOrientation','portrait',...
    'position',[400 100   800  400]);
names = {'vp1';'vp2';'vp3'};
symb = {'o';'d';'p'};
lines = {'-';'-';'-'};
tstep = datenum(0,0,0,0,60,0);
hight = datenum(2015,03,09,04,33,00);
sp = zeros(4,1);
cl = flipud([0.1 0.1 0.1;0.5 0.5 0.5;0.7 0.7 0.7]);
pl = zeros(3,1);
space = [300 250 400];
for i = 1:3
    time = data.(names{i}).time;
    lvl = data.(names{i}).ble.*100;
%     if i == 1
%         lvl = smooth(lvl,50);
%     end
    plot(time,zeros(length(time),1),'-k'),hold on
    plot(hight*ones(10,1),linspace(-10,10,10),'--k',...
        'linewidth',1.5)
    plot(time,lvl,lines{i},...
        'color',cl(i,:),'linewidth',1.5)
    lc = 1:space(i):length(time);
    pl(i) = plot(time(lc),lvl(lc),symb{i},...
        'LineWidth',1.5,...
        'markeredgecolor','k',...
        'markerfacecolor',cl(i,:),...
        'markersize',8);
    if i == 3
        datetick('x','HH:MM','keepticks','keeplimits')
        set(gca,'xlim',[time(1) time(end)],...
        'xtick',time(1):tstep:time(end),...
        'ylim',[-5 2])
    end
end
leg = legend(pl,{'x = -10cm';'x = 10cm';'x = 20cm'});
set(leg,'position',[0.85 0.3 0.05 0.05])
ylabel('Bed Level Change (cm)')
xlabel('Time on 09/03/15')
% prettyfigures('text',13,'labels',14,'box',1)
% export_fig([figdir 'HTABottomTrackTS'],'-pdf')
break
%Load bed change data
f2 = figure(2);
set(f2,'PaperOrientation','portrait',...
    'position',[400 100   800   400]);
dirc = 'd:\Projects\Mekong_W2015\Images\Quadrats\Q2\Reconstructions\';
qdat = 'D:\Projects\Mekong_W2015\Images\Quadrats\Q2\Reconstructions\2A\Analysis\';
fname = {'twoatwobdiff.txt','twobtwocdiff.txt'};
qname = {'Q2Astats_fixed.mat','Q2Bstats_fixed.mat'};
vpxy = [0.8033 0.245;0.5997 0.245;0.4978 0.245];
rotx = [0,0];
roty = [0,0];
rotz = [0,0];
cc = brewermap(100,'RdYlBu');
colormap(cc);
for i = 1:2
    [z,r] = arcgridread([dirc fname{i}]);
    ccscale = [7.29 6.657];
    Z = (z./ccscale(i)).*100; %*100 converts to cm
%     Z = rot90(Z,-1); %rotate Z so it is in the same view as the orthophotos
    %create unit vector axis
    [m,n] = size(Z);
    x = linspace(0,1,n);X = repmat(x,m,1);
    y = linspace(1,0,m);Y = repmat(y,n,1)';
    %load quadrat data
    load([qdat 'Q2Astats_fixed.mat'])
    
    %contour using contourf
    sp(i) = subplot(1,2,i);
    [c,h] = contourf(X,Y,Z);
    set(h,'EdgeColor','none')
    hold on
    %calculate radii of pneumatophores at level 1 of the reconstruction
    nn = length(DATA.layer1.Rcenter);
    for j = 1:nn %loop through rows for the Radius and Center points of each layer
        cx = DATA.layer1.Rcenter(j,1);
        cy = DATA.layer1.Rcenter(j,2);
        %transformation matrix, rotate 90deg CCW
        R = [0 -1;1 0];
        rot = R*[cx;cy];
        cx = abs(rot(1));cy = rot(2);
        cr = mean(DATA.layer1.Rradii(j,:));
        ang=0:0.01:2*pi;
        xp=cr.*cos(ang);
        yp=cr.*sin(ang);
        %plot circles with the same radius as the vegetation
        H = patch(cx+xp,cy+yp,1);
        set(H,'FaceColor',[0 0 0],'EdgeColor','none')
    end
    for j = 1:3
        plot(vpxy(j,1),vpxy(j,2),symb{j},...
            'Color','k',...
            'MarkerSize',10,...
            'MarkerFaceColor',cl(j,:),...
            'LineWidth',1.5)
    end
    caxis([-3 3])
    grid on
    set(gca,'gridlinestyle',':')
end
%labels
cb = colorbar;
set(cb,'ytick',-3:1.5:3)
ylabel(cb,'Bed Level Change (cm)')
ylabel(sp(1),'Along Shore (m)')
xlabel(sp(1),'Cross Shore (m)')
xlabel(sp(2),'Cross Shore (m)')
%positioning
%positioning
set(sp(1),'position',[0.1 0.14 0.36 0.8],...
    'xtick',0:0.25:1,'ytick',0:0.25:1)
set(sp(2),'position',[0.5 0.14 0.36 0.8],...
    'yticklabel',[],'xtick',0:0.25:1,...
    'ytick',0:0.25:1)
set(cb,'position',[0.89 0.14 0.02 0.8])
prettyfigures('text',13,'labels',14,'box',1)
% export_fig([figdir 'HTABedChange'],'-pdf')