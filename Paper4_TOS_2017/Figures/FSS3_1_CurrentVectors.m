%Figure out the current rotations for the FSS_1 (2015) experiment

close all
clear

pw = 'D:\Projects\Mekong_W2015\Data\Aquadopp\FSS\';
files = dir([pw '*9March2015.mat']);files = {files.name};
files{3} = [];files{4} = [];
names = {'5116';'5117';'V5109';'VC01';'VC101'};
symb = {'o';'^';'d';'p';'h'};
scale = 0.01;
start = datenum(2015,03,07,14,30,00);
stop = datenum(2015,03,07,17,38,00);
count = 1;
for ff = 1:length(files)
    if isempty(files{ff})
        continue
    else
        load([pw files{ff}])
        disp(['Loading ' files{ff}])
        
        ind = find(aqdp.datenum >= start & aqdp.datenum <= stop);
        if ff < 3
            u = cmgbridge(nanmean(aqdp.u(ind,:),2),100,1000,10000);
            v = cmgbridge(nanmean(aqdp.v(ind,:),2),100,1000,10000);
        else
            u = cmgbridge(nanmean(aqdp.u(ind,:),2),100,1000,10000);
            v = cmgbridge(nanmean(aqdp.v(ind,:),2),100,1000,10000);
        end
        u(isnan(u)) = 0;v(isnan(v)) = 0;
        lat = aqdp.metadata.lat;
        lon = aqdp.metadata.lon;
        time = aqdp.datenum(ind);
        pres = cmgbridge(aqdp.pressure(ind),100,1000,10000);
        fs = 8;
        intv = 7.5;
        avt  = fs*intv*60;
        id = [1 avt:avt:length(u)];
        
        figure(1)
        if ff < 3
            c = jet(length(id));
            for i = 1:length(id)-1
                quiver(lon,lat,mean(u(id(i):id(i+1))),mean(v(id(i):id(i+1))),scale,...
                    'Color',c(i,:));hold on
            end
        else
            c = jet(length(u));
            for i = 1:length(u)
                quiver(lon,lat,u(i),v(i),scale,...
                    'Color',c(i,:));hold on
            end
        end
        text(lon,lat,names{count},'FontSize',10)
        
        figure(2)
        if ff < 3
            c = jet(length(id));
            for i = 1:length(id)-1
                pp(count) = plot(time(id(i)),mean(pres(id(i):id(i+1))),...
                    'Color','k','MarkerFaceColor',c(i,:),...
                    'Marker',symb{count});hold on
            end
        else
            c = jet(length(u));
            for i = 1:length(u)
                pp(count) = plot(time(i),pres(i),...
                    'Color','k','MarkerFaceColor',c(i,:),...
                    'Marker',symb{count});hold on
            end
        end
        datetickzoom('x','HH:MM:SS','keepticks','keeplimits')
        
        figure(3)
        n = length(u);
        pd(count) = plot(lon,lat,'Marker',symb{count},...
            'MarkerSize',8,'LineWidth',1,'Color','k'); hold on
        if ff < 3
            c = jet(length(id));
            for i = 1:length(id)-1
                umean(i) = mean(u(id(i):id(i+1)));
                vmean(i) = mean(v(id(i):id(i+1)));
                posX = cumsum([lon ; umean(:).*scale]);
                posY = cumsum([lat ; vmean(:).*scale]);
            end
            for i = 1:length(id)-1
                line(posX(i:i+1),posY(i:i+1),...
                    'LineWidth',1.5,...
                    'Color',c(i,:));
            end
        else
            u(isnan(u)) = 0;v(isnan(v)) = 0;
            c = jet(n);
            posX = cumsum([lon ; u(:).*scale]);
            posY = cumsum([lat ; v(:).*scale]);
            for i = 1:n-1
                line(posX(i:i+1),posY(i:i+1),...
                    'LineWidth',1.5,...
                    'Color',c(i,:));
            end
        end
        count = count+1;
    end
end
clear aqdp
%Load the Vector
pw = 'D:\Projects\Mekong_W2015\Data\Vector\FSS\';
files = dir([pw '*070315.mat']);files = {files.name};

for ff = 1:length(files)
    load([pw files{ff}])
    disp(['Loading ' files{ff}])
    
    ind = find(ADV.datetime >= start & ADV.datetime <= stop);
    u = cmgbridge(nanmean(ADV.U(ind),2),100,1000,10000);
    v = cmgbridge(nanmean(ADV.V(ind),2),100,1000,10000);
    u(isnan(u)) = 0;v(isnan(v)) = 0;
    lat = ADV.Metadata.inst_lat;
    lon = ADV.Metadata.inst_lon;
    time = ADV.datetime(ind);
    pres = cmgbridge(ADV.Pres(ind),100,1000,10000);
    fs = 32;
    intv = 8;
    avt  = fs*intv*60;
    id = [1 avt:avt:length(u)];
    
    f1 = figure(1);
    c = jet(length(id));
    for i = 1:length(id)-1
        quiver(lon,lat,mean(u(id(i):id(i+1))),mean(v(id(i):id(i+1))),scale,...
            'Color',c(i,:));hold on
    end
    
    text(lon,lat,names{count},'FontSize',10)
    
    f2 = figure(2);
    c = jet(length(id));
    for i = 1:length(id)-1
        pp(count) = plot(time(id(i)),mean(pres(id(i):id(i+1))),...
            'Color','k','MarkerFaceColor',c(i,:),...
            'Marker',symb{count});hold on
    end
    
    %     datetickzoom('x','HH:MM:SS','keepticks','keeplimits')
    
    f3 = figure(3);
    n = length(u);
    pd(count) = plot(lon,lat,'Marker',symb{count},...
        'MarkerSize',8,'LineWidth',1,'Color','k'); hold on
    c = jet(length(id));
    for i = 1:length(id)-1
        umean(i) = mean(u(id(i):id(i+1)));
        vmean(i) = mean(v(id(i):id(i+1)));
        posX = cumsum([lon ; umean(:).*scale]);
        posY = cumsum([lat ; vmean(:).*scale]);
    end
    for i = 1:length(id)-1
        line(posX(i:i+1),posY(i:i+1),...
            'LineWidth',1.5,...
            'Color',c(i,:));
    end
    count = count+1;
end
savefigdir = 'd:\Projects\Mekong_W2015\Figures\Environment\';
figure(1)
xlabel('lon')
ylabel('lat')
title('Current Vectors, FSS3-1')
export_fig([savefigdir 'FSS3_1_allCurvec'],'-pdf','-nocrop')
figure(2)
xlabel('Time on 07/03/2015')
ylabel('dbar')
title('Pressure Records, FSS3-1 Multiple Instruments')
legend(pp,names)
export_fig([savefigdir 'FSS3_1_allPres'],'-pdf','-nocrop')
figure(3)
xlabel('lon')
ylabel('lat')
title('Progressive Vectors, FSS3-1 Multiple Instruments')
legend(pd,names)
export_fig([savefigdir 'FSS3_1_allPVDs'],'-pdf','-nocrop')




