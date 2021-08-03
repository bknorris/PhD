%This script plots cospectra from all 3 days of the HTA and single day of
%the VTA, and calculates the phase between the two signals (x-shore
%velocity & turbulence).

clear
close all

%%%Load Data\Data prep for plotting
datdir = 'd:\Projects\Mekong_W2015\DataAnalysis\Paper1\Eps_Vel_Spectra\';
load('D:\Projects\Mekong_W2015\Data\Vectrino\VPTmatrix.mat')
load([datdir 'VectorVels.mat'])
velfiles = dir([datdir '*30min.mat']);
velfiles = {velfiles.name};
tkefiles = dir([datdir 'Turbulence\' '*30minTKE.mat']);
tkefiles = {tkefiles.name};
PCA = 1;
plotfigs = 1;
timeseries = 1;
savefigdir = 'd:\Projects\Mekong_W2015\Figures\Paper1\WaveVelsTurbulence\';
data = struct();
bins = [5 5 5;15 15 15;15 15 15;5 15 15];
for f = 1:4
    disp(['Loading ' velfiles{f}])
    disp(['Loading ' tkefiles{f}])
    load([datdir velfiles{f}])
    load([datdir 'Turbulence\' tkefiles{f}])
    vname = velfiles{f}(1:5);
    if strcmp(vname(1:3),'HTA')
        name = [vname(1:3) num2str(f)];
    else
        name = regexprep(vname,'_','');
    end
    fn = fieldnames(dat);
    tsl = zeros(3,1);
    %find the shortest t-s; sometimes the instrument datasets are not
    %the same length
    for inst = 1:3
%         T = Tmat.(fn{inst});
%         b1 = dat.(fn{inst}).beam1;b2 = dat.(fn{inst}).beam2;
%         b3 = dat.(fn{inst}).beam3;b4 = dat.(fn{inst}).beam4;
%         [x,y,z1,z2] = VPTransform(b1,b2,b3,b4,T,'bx');
        tsl(inst) = length(dat.(fn{inst}).x);
%         clear b1 b2 b3 b4
%         dat.(fn{inst}).x = x;dat.(fn{inst}).y = y;
%         dat.(fn{inst}).z1 = z1;dat.(fn{inst}).z2 = z2;
    end
    [tsl,id] = min(tsl);
    start = dat.(fn{id}).time(1);stop = dat.(fn{id}).time(end);
    %timebases are off by several milliseconds between instruments.
    %Create new timebase
    time = linspace(start,stop,tsl);
    %use PCA to calculate the mean wave/current direction from
    %ancillary instrument
    if PCA == 1
        u = vector.(name).u;v = vector.(name).v;
        vid = find(vector.(name).time >= start & vector.(name).time <= stop);
        u = cmgbridge(u(vid),100,100,1000);v = cmgbridge(v(vid),100,100,1000);
        cmp = pca(u,v,[],0);heading = cmp.mdir;
        %VPs point to N or 20deg (HTA/VTA)
        if strcmp(name(1:3),'HTA')
            heading = 360-abs(0-heading);
        else
            heading = 360-abs(20-heading);
        end
    else
        heading = 20;
    end
    
    for inst = 1:3 %loop through vectrinos
        bin = bins(f,inst);
        x = dat.(fn{inst}).x(:,bin);y = dat.(fn{inst}).y(:,bin);
                
        %rotate to cross-shore and along-shore
        rot = heading*pi/180;
        x = x.*(ones(size(x))*cos(rot)) + ...
            y.*(ones(size(y))*sin(rot));
        y = -y.*(ones(size(y))*sin(rot)) + ...
            x.*(ones(size(x))*cos(rot));
        
        %average per-beam turbulence statistics into a single time-series
        eps = (Stat.(fn{inst}).beam1.E(bin,:)+Stat.(fn{inst}).beam2.E(bin,:)...
            +Stat.(fn{inst}).beam3.E(bin,:)+Stat.(fn{inst}).beam4.E(bin,:))./4;
        eps = cmgbridge(eps,10,100,100);eps(isnan(eps)) = 0;
        time1 = Stat.(fn{inst}).time;
        
        %interpolate turbulence (@ 1 Hz) up to 50 Hz
        time2 = dat.(fn{inst}).time;
        eps = spline(time1,eps,time2);
        eps(eps < 0) = 0; %set any negative values to zero
        
        %crop to shortest t-s
        x = x(1:tsl);y = y(1:tsl);eps = eps(1:tsl);
        
        %detrend
        x = detrend(x);y = detrend(y);
        
        %compute cospectra
        fs  = 50; %Hz
        step = 300; %seconds
        avt = fs*step; %# samples per window
        window = avt*0.2;
        nfft = window*0.5; %equal to 30 second subsamples
        idx = [1 avt:avt:tsl];
        n = length(idx)-1;
        Sxy = zeros(nfft/2+1,length(idx)-1);
        Cxy = zeros(nfft/2+1,length(idx)-1);
        for i = 1:n
            xx = x(idx(i):idx(i+1));
            epsx = eps(idx(i):idx(i+1));
            [Sxy(:,i),F] = cpsd(xx,epsx,hanning(window),window*0.7,nfft,fs);
            [Cxy(:,i),~] = mscohere(xx,epsx,hanning(window),window*0.7,nfft,fs);
        end
        %magnitudes in dB
        samp = Sxy.*conj(Sxy);
        camp = Cxy.*conj(Cxy);
        smag = 10*log10(abs(Sxy));
        phase = angle(Sxy);
        
        %save data into structure
        data.(name).(fn{inst}).time = time;
        data.(name).(fn{inst}).x = x;
        data.(name).(fn{inst}).y = y;
        data.(name).(fn{inst}).eps = eps;
        data.(name).(fn{inst}).F = F;
        data.(name).(fn{inst}).Sxy = Sxy;
        data.(name).(fn{inst}).Cxy = Cxy;
        data.(name).(fn{inst}).Camp = camp;
        data.(name).(fn{inst}).Samp = samp;
        data.(name).(fn{inst}).phase = phase;
    end
    clear dat Stat x y eps
end
if plotfigs
    %%%Plot Routine
    f1 = figure(1);
    set(f1,'PaperOrientation','portrait',...
        'position',[400 200   1200   600]);
    set(gcf,'color','w','PaperPositionMode','auto')
    sp = zeros(3,1);p = zeros(3,1);
    dn = fieldnames(data);
    symb = {'^';'o';'s'};
    c = [0.6 0.6 0.6;0.4 0.4 0.4;0 0 0];
    row = zeros(4,3);col = zeros(4,3);
    for i = 1:4
        sp(i) = subplot(1,4,i);
        for j = 1:3
            F = data.(dn{i}).(fn{j}).F;
            S = data.(dn{i}).(fn{j}).Samp;
            phase = data.(dn{i}).(fn{j}).phase;
            F(1) = [];S(1,:) = [];phase(1,:) = []; %remove DC offset
            %find maximum value; plot as an example
            [~,id] = max(max(S));[~,id2] = max(S(:,id));
            p(j) = plot(F,S(:,id),'color',c(j,:),...
                'LineWidth',1.5);hold on
            plot(F(id2),S(id2,id),symb{j},...
                'markeredgecolor','k',...
                'markerfacecolor',c(j,:))
            text(F(id2)+0.05,S(id2,id),['(' num2str(F(id2),2) ', ' num2str(S(id2,id),3) ')'])
            text(0.15,5E-3,dn{i})
            %calculate phase
            theta = rad2deg(phase(id2,id));
            f = F(id2);
            deltaT = theta/(360*f);
            disp([dn{i} ' ' fn{j} ' time lag: ' num2str(deltaT) ' seconds'])
            row(i,j) = id2;col(i,j) = id;
        end
        
    end
    %global adjustments
    legend(p,fn)
    labeltext = {'(Frequency, Peak)'};
    axes(sp(1))
    text(0.15,1E-3,labeltext)
    h1 = ylabel(sp(1),'Power');
    [~,h2] = suplabel('Frequency (Hz)','x',[.1 .08 .84 .84]);
    [~,t1] = suplabel('Cross Spectra - All Experiments','t');
    set([h1 h2 t1],'FontSize',14,'FontName','Arial')
    set([sp(2) sp(3) sp(4)],'YTickLabel',[])
    set(sp,'XLim',[0 1.2],'yscale','log','ylim',[1E-15 1E-2])
    if PCA
        pc = 'PCA';
    else
        pc = '20deg';
    end
    export_fig([savefigdir 'Cospectra\' 'CPSD_allExp_' pc],'-png')
    
    f2 = figure(2);
    set(f2,'PaperOrientation','portrait',...
        'position',[400 200   1200   600]);
    set(gcf,'color','w','PaperPositionMode','auto')
    sp = zeros(3,1);p = zeros(3,1);
    dn = fieldnames(data);
    symb = {'^';'o';'s'};
    c = [0.6 0.6 0.6;0.4 0.4 0.4;0 0 0];
    for i = 1:4
        sp(i) = subplot(1,4,i);
        for j = 1:3
            F = data.(dn{i}).(fn{j}).F;
            C = data.(dn{i}).(fn{j}).Camp;
            phase = data.(dn{i}).(fn{j}).phase;
            F(1) = [];C(1,:) = [];
            range = find(F > 0 & F < 1.2);
            F = F(range);C = C(range,:);
            %         [~,id] = max(max(C));[~,id2] = max(C(:,id));
            id2 = row(i,j);id = col(i,j);
            p(j) = plot(F,C(:,id),'color',c(j,:),...
                'LineWidth',1.5);hold on
            plot(F(id2),C(id2,id),symb{j},...
                'markeredgecolor','k',...
                'markerfacecolor',c(j,:))
            text(F(id2)+0.05,C(id2,id),['(' num2str(F(id2),2) ', ' num2str(C(id2,id),3) ')'])
            text(0.15,0.95,dn{i})
        end
    end
    %global adjustments
    legend(p,fn)
    labeltext = {'(Frequency, MSC)';'D.O.F. = 20';...
        '95% Significance: 0.146';'99.9% Significance: 0.305'};
    axes(sp(1))
    text(0.15,0.85,labeltext)
    prettyfigures('font','arial','text',14,'labels',10,'box',1,'grid',0)
    h1 = ylabel(sp(1),'\gamma^2');
    [~,h2] = suplabel('Frequency (Hz)','x',[.1 .08 .84 .84]);
    [~,t1] = suplabel('MSC - All Experiments','t');
    set([h1 h2 t1],'FontSize',14,'FontName','Arial')
    set([sp(2) sp(3) sp(4)],'YTickLabel',[])
    set(sp,'XLim',[0 1.2],'ylim',[0 1])
    if PCA
        pc = 'PCA';
    else
        pc = '20deg';
    end
    export_fig([savefigdir 'Cospectra\' 'MSC_allExp_' pc],'-png')
    
    if timeseries
        close('all')
        %%%Plot snapshot of velocity & turbulence & phase lag of turbulence
        %%%behind velocity
        times = [datenum(2015,03,07,15,45,00) datenum(2015,03,07,15,47,00);...
            datenum(2015,03,08,15,38,52) datenum(2015,03,08,15,40,52);...
            datenum(2015,03,10,16,19,12) datenum(2015,03,10,16,22,12);...
            datenum(2015,03,14,08,55,40) datenum(2015,03,14,08,57,40)];
        ylims = [0 2.5E-2;0 3E-2;0 1.5E-2;0 4E-3];
        for i = 1:4
            for j = 1:3
                x = data.(dn{i}).(fn{j}).x;
                y = data.(dn{i}).(fn{j}).y;
                eps = data.(dn{i}).(fn{j}).eps;
                time = data.(dn{i}).(fn{j}).time;
                p = zeros(2,1);
                if strcmp(fn{j},'vpro1') && strcmp(dn{i}(1:3),'HTA')
                    id = find(x<0);x2 = NaN(length(x),1);x2(id) = x(id);
                    id = find(y<0);y2 = NaN(length(x),1);y2(id) = y(id);
                else
                    id = find(x>0);x2 = NaN(length(x),1);x2(id) = x(id);
                    id = find(y>0);y2 = NaN(length(x),1);y2(id) = y(id);
                end
                
                f1 = figure;
                set(f1,'PaperOrientation','portrait',...
                    'position',[400 200   800   500]);
                set(gcf,'color','w','paperpositionmode','auto')
                sp(1) = subplot(211);
                p(1) = plot(time,x,'-','color','k',...
                    'LineWidth',1.5);hold on
                plot(time,x2,'-','color','g',...
                    'linewidth',1.5)
                p(2) = plot(time,y,'-.','color',[0.5 0.5 0.5],...
                    'linewidth',1.5);
                plot(time,y2,'-.','color','r',...
                    'linewidth',1.5)
                plot(time,zeros(length(x),1),'Color',[0.5 0.5 0.5])
                set(gca,'xlim',[times(i,1) times(i,2)],'Xticklabel',[],...
                    'ylim',[-0.8 0.8],'ytick',-0.8:0.4:0.8)
                title(['Unjusted Velocities - ' upper(fn{j}) ', ' dn{i}])
                ylabel('Velocity (m/s)')
                legend(p,{'x-shore';'along-shore'},...
                    'position',[0.87 0.83 0.05 0.05])
                sp(2) = subplot(212);
                p(1) = plot(time,eps,'k','linewidth',1.5); hold on
                set(gca,'xlim',[times(i,1) times(i,2)],...
                    'ylim',[ylims(i,1) ylims(i,2)])
                linkaxes(sp,'x')
                datetickzoom('x','HH:MM:SS','keepticks','keeplimits')
                xlabel(['Time on ' datestr(time(1),'dd-mm-yyyy')])
                ylabel('\epsilon (Wkg^-^1')
                title(['Adjusted Turbulence - ' upper(fn{j}) ', ' dn{i}])
                tname = datestr(times(i,1),'HHMMSS');
                export_fig([savefigdir 'PCAvels\' dn{i} '_' fn{j} '_' tname],'-png')
            end
        end
    end
end




