clear
datdir = 'c:\Users\bkn5\Projects\Mekong_W2015\DataAnalysis\Spectra\Paper1\QCd\';
figdir = 'C:\Users\bkn5\Projects\Mekong_W2015\Figures\Velocity\';
load([datdir 'HTAday1Vels']);
load('C:\Users\bkn5\Projects\Mekong_W2015\DataAnalysis\VPs\Paper1\HTAday1BdTrack.mat')
%time average velocities:
start = datenum(2015,03,07,14,15,00);stop = datenum(2015,03,07,17,10,00);tstep = datenum(0,0,0,0,30,0);
fn = fieldnames(HTA);bfn = fieldnames(dat);
hab = [0.063 0.061 0.062]*100;
vname = {'VP1';'VP2';'VP3'};
count = 1;
for i = 2:4
    %velocities
    idx = find(HTA.(fn{i}).time >= start & HTA.(fn{i}).time <= stop);
    h = hab(count); %inst hab (cm)
    ph = 43.5; %pneumatophore height (cm)
    
    u = HTA.(fn{i}).y(idx,1:end); %along-shore
    v = HTA.(fn{i}).x(idx,1:end); %cross-shore
    w = HTA.(fn{i}).z(idx,1:end);
    t = HTA.(fn{i}).time(idx,:);
    
    %bottom distance
    idx = find(dat.(bfn{count}).bchtm >= start & dat.(bfn{count}).bchtm  <= stop);
    bdist = dat.(bfn{count}).bdist(idx,:);
    bchtm = dat.(bfn{count}).bchtm(idx,:);
    bdmax = runningmax(bdist,100);
    bdmax = my_running_median(bdmax,500); %despike
    bdmax = smooth(bdmax,100,'sgolay');
    %window to smooth the line further
    win = 400; 
    step = 200;
    ind = [1 step:step:length(bdmax)];
    bdavs = zeros(length(ind),1);
    bdavt = zeros(length(ind),1);
    for ii = 1:length(ind)
        if abs(length(bdmax)-ind(ii)) < win
            continue
        else
            bwin = bdmax(ind(ii):ind(ii)+win);
            bdavs(ii,:) = nanmean(bwin);
            bdavt(ii,:) = bchtm(ind(ii));
        end
    end
    bdavs(bdavs == 0) = []; %remove trailing zeros
    bdavt(bdavt == 0) = [];
    bdavs = h-bdavs.*100;
    
    rb = HTA.(fn{i}).rb;rbcm = h-rb(1:end)*100;
    sr = 50;
    intv = 2; %minutes
    avt = intv*sr*60;
    ind = [1 avt:avt:length(u)];
    
    f1 = figure(i);
    set(f1,'PaperOrientation','portrait',...
        'position',[200 200   1000   600]);
    set(gcf,'color','w','PaperPositionMode','auto')
    p(1) = subplot(131);
    for ii = 1:length(ind)-1
        Uz = u(ind(ii):ind(ii+1),:);
        U(ii,:) = nanmean(Uz);
        T(ii,:) = t(ind(ii));
    end
    % c = jet(length(ind));
    imagesc(T,rbcm,U'*100), hold on
    plot(bdavt,bdavs,'-k')
    datetickzoom('x','HH:MM','keepticks','keeplimits')
    set(gca,'Ydir','normal','YGrid','on','xlim',[T(1) T(end)],...
        'xtick',T(1):tstep:T(end),'Ylim',[-1 2.5])
    caxis([-2.5 2.5])
    xlabel('Velocity (cm/s)')
    ylabel('Height Above Bed (cm)')
    title('Along-Shore')
    
    p(2) = subplot(132);
    for ii = 1:length(ind)-1
        Vz = v(ind(ii):ind(ii+1),:);
        V(ii,:) = nanmean(Vz);
    end
    imagesc(T,rbcm,V'*100), hold on
    plot(bdavt,bdavs,'-k')
    datetickzoom('x','HH:MM','keepticks','keeplimits')
    set(gca,'Ydir','normal','YGrid','on','xlim',[T(1) T(end)],...
        'xtick',T(1):tstep:T(end),'Ylim',[-1 2.5],'YTickLabel',[])
    caxis([-2.5 2.5])
    xlabel('Velocity (cm/s)')
    % ylabel('Height Above Bed (cm)')
    title('Cross-Shore')
    
    p(3) = subplot(133);
    for ii = 1:length(ind)-1
        Wz = w(ind(ii):ind(ii+1),:);
        W(ii,:) = nanmean(Wz);
    end
    imagesc(T,rbcm,W'*100), hold on
    plot(bdavt,bdavs,'-k')    
    datetickzoom('x','HH:MM','keepticks','keeplimits')
    set(gca,'Ydir','normal','YGrid','on','xlim',[T(1) T(end)],...
        'xtick',T(1):tstep:T(end),'Ylim',[-1 2.5],'YTickLabel',[])
    caxis([-2.5 2.5])
    xlabel('Velocity (cm/s)')
    % ylabel('Height Above Bed (cm)')
    title('Vertical')
    
    set(p(1),'position',[0.1 0.22 0.25 0.7])
    set(p(2),'position',[0.38 0.22 0.25 0.7])
    set(p(3),'position',[0.66 0.22 0.25 0.7])
    cb = colorbar('south');
    set(cb,'position',[0.25 0.02, 0.54 0.05]),xlabel(cb,'cm/s')
    suptitle([vname{count} ' : ' datestr(start,'dd/mm HH:MM')])
    
    %save figs
    export_fig([figdir 'HTAday1' vname{count} 'depVels'],'-jpeg')
    count = count+1;
end
