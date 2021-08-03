%Plot bed movement during the HTA experiments (A re-do of the original
%'scour.m' plot
clear
close all
figdir = 'd:\Projects\Mekong_W2015\Figures\Paper1\';
vph = [0.062 0.063 0.061 0.242 0.271 0.240];
datdir = 'D:\Projects\Mekong_W2015\Data\Vectrino\7March2015\';
fname = dir([datdir '*070315.mat']);
inst = cell(6,1);
for i = 1:3
    inst{i} = [datdir char({fname(i).name})];
end
datdir = 'D:\Projects\Mekong_W2015\Data\Vectrino\8March2015\';
fname = dir([datdir '*080315.mat']);
ct = 4:6;
for i = 1:3
    inst{ct(i)} = [datdir char({fname(i).name})];
end

datdir = 'd:\Projects\Mekong_W2015\DataAnalysis\Paper1\';
if exist([datdir 'BDtrace_HTA1_2.mat'],'file') > 0
    disp('Loading BDtrace_HTA1_2.mat')
    load([datdir 'BDtrace_HTA1_2.mat'])
else
    for i = 1:6
        if i < 4
            date = 'day1';
            start = datenum(2015,03,07,13,36,00);stop = datenum(2015,03,07,17,10,00);
        else
            date = 'day2';
            start = datenum(2015,03,08,14,15,00);stop = datenum(2015,03,08,18,30,00);
        end
        vname = lower(inst{i}(51:53));
        
        disp(['Loading ' inst{i}])
        load(inst{i})
        bd = VPRO.Data.BottomCheck_BottomDistance;
        gmt2ict = datenum(0,0,0,1,0,0)*7;
        time = VPRO.Data.BottomCheck_HostTimeMatlab+gmt2ict;
        id = find(time >= start & time <= stop);
        bd = bd(id);time = time(id);
        
        %%%Smooth bottom distance
        bdmid = my_running_median(bd,100);
        %window to smooth for plotting
        win = 1000;
        step = 100;
        idx = [1 step:step:length(bdmid)];
        bdmax = zeros(length(idx),1);
        time2 = zeros(length(idx),1);
        for ii = 1:length(idx)-1
            if abs(length(bdmid)-idx(ii)) < win
                continue
            else
                bwin = bdmid(idx(ii):idx(ii)+win);
                bdmax(ii,:) = max(bwin);
                time2(ii,:) = time(idx(ii));
            end
        end
        bdmax(bdmax == 0) = []; %remove trailing zeros
        time2(time2 == 0) = [];
        %     id2 = find(bdmax < mean(bdmax)-2.5*std(bdmax));
        %     if ~isempty(id2)
        %         if id2(1) < 21
        %             rm = min(id2):max(id2)+20;
        %         elseif max(id2) > length(bdmax)-21
        %             rm = min(id2)-20:max(id2);
        %         else
        %             rm = min(id2)-20:max(id2)+20;
        %         end
        %         bdmax(rm) = NaN;
        %         bad=isnan(bdmax);
        %         gd=find(~bad);
        %         bad([1:(min(gd)-1) (max(gd)+1):end])=0;
        %         bdmax(bad)=interp1(gd,bdmax(gd),find(bad),'pchip');
        %     end
        if i < 4
            bdavs = smooth(bdmax,25);
        else
            bdavs = smooth(bdmax,100);
        end
        data.(date).(vname).time = time2;
        data.(date).(vname).bd = bdavs;
        data.(date).(vname).lvl = vph(i)-bdavs;
        clear VPRO
    end
    save([datdir 'BDtrace_HTA1_2'],'data','-v7.3')
end
%%%Plot Routine%%%
f1 = figure(1);
set(f1,'PaperOrientation','portrait',...
    'position',[400 100   1200   800]);
names = {'vp1';'vp2';'vp3'};
symb = {'o';'d';'p'};
lines = {'-';'--';'-.'};
tstep = datenum(0,0,0,0,30,0);
sp = zeros(4,1);
cl = flipud([0.1 0.1 0.1;0.5 0.5 0.5;0.7 0.7 0.7]);
pl = zeros(3,1);
for i = 1:3
    sp(1) = subplot(221);
    time = data.day1.(names{i}).time;
    lvl = data.day1.(names{i}).lvl.*100;
    plot(time,zeros(length(time),1),'-k'),hold on
    plot(time,lvl,lines{i},...
        'color',cl(i,:),'linewidth',1.5)
    lc = 1:200:length(time);
    plot(time(lc),lvl(lc),symb{i},...
        'LineWidth',1.5,...
        'markeredgecolor','k',...
        'markerfacecolor',cl(i,:),...
        'markersize',8)
    set(sp(1),'xlim',[time(1) time(end)],...
        'xtick',time(1):tstep:time(end),...
        'ylim',[-5 2])
    datetick('x','HH:MM','keepticks','keeplimits')
    sp(2) = subplot(222);
    time = data.day2.(names{i}).time;
    lvl = data.day2.(names{i}).lvl.*100;
    if i == 2
        lvl(2:23) = NaN;
        bad=isnan(lvl);
        gd=find(~bad);
        bad([1:(min(gd)-1) (max(gd)+1):end])=0;
        lvl(bad)=interp1(gd,lvl(gd),find(bad),'pchip');
    end
    plot(time,zeros(length(time),1),'-k'),hold on
    plot(time,lvl,lines{i},...
        'color',cl(i,:),'linewidth',1.5)
    pl(i) = plot(time(lc),lvl(lc),symb{i},...
        'LineWidth',1.5,...
        'markeredgecolor','k',...
        'markerfacecolor',cl(i,:),...
        'markersize',8);
    set(sp(2),'xlim',[time(1) time(end)],...
        'xtick',time(1):tstep:time(end),...
        'ylim',[-10 2])
    datetick('x','HH:MM','keepticks','keeplimits')
end
leg = legend(pl,{'x = -10cm';'x = 10cm';'x = 20cm'});
ylabel(sp(1),'Bed Level Change (cm)')
xlabel(sp(1),'Time on 07/03/15')
xlabel(sp(2),'Time on 08/03/15')

%Load bed change data
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
    sp(i+2) = subplot(2,2,i+2);
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
ylabel(sp(3),'Along Shore (m)')
xlabel(sp(3),'Cross Shore (m)')
xlabel(sp(4),'Cross Shore (m)')
title(sp(1),'HTA Day 1')
title(sp(2),'HTA Day 2')
%positioning
set(sp(1),'position',[0.1 0.56 0.34 0.38])
set(sp(2),'position',[0.5 0.56 0.34 0.38])
set(sp(3),'position',[0.1 0.1 0.34 0.38])
set(sp(4),'position',[0.5 0.1 0.34 0.38])
set(leg,'position',[0.8 0.6 0.05 0.05])
set(cb,'position',[0.855 0.1 0.02 0.38])
prettyfigures('text',13,'labels',14,'box',1)

export_fig([figdir 'HTAbottomTrack'],'-pdf')