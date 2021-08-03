%Calculate the slope of velocity spectra from VecPros during the VTA and
%VTA. Calculate spectra at set heights: Bins 1,5,10,15,20,25,35 for each
%VecPro during the experiments. Capture the slope of the spectra from 1.2Hz
%to 10Hz, and save it to a variable.
clear
close all
ploton = 0;
datdir = 'c:\Users\bkn5\Projects\Mekong_W2015\DataAnalysis\Paper1\';
Spec = struct();
heading = 100; %rotate vecpros pointing 20 deg CCW + (360-280 (xshore angle))
%define instrument heights from deployment logs
H.vpro1 = 0.07;H.vpro2 = 0.416;H.vpro3 = 0.806;
%Check for the existence of the datafile to conditionally run the program
ifspec = exist([datdir 'VTASlopeSpectra.mat'],'file');
disp('Loading VTAvelocities.mat')
load([datdir 'VTAvelocities.mat'])
S = struct();
%loop through the instruments
for ii = 2:4
    fn = fieldnames(VTA);
    disp(['Running spectral analysis on ' fn{ii}])
    hab = H.(fn{ii});
    fs = VTA.(fn{ii}).sr;
    
    %window the time series. For now, use a 30sec window and 10 second step
    win = 30; %second
    step = 10; %second
    nsamp = win*fs;
    ind = 1:step*fs:length(VTA.(fn{ii}).x);
    for j = 1:length(ind)
        if abs(length(VTA.(fn{ii}).x)-ind(j)) < nsamp %skip the last few indexes approaching the end of the t-s
            continue
        else
            idx = ind(j):ind(j)+nsamp; %steps include number of samples equal to the window size * the sampling freq
            x = VTA.(fn{ii}).x(idx,:);
            y = VTA.(fn{ii}).y(idx,:);
            %rotate to cross-shore:
            rot = heading*pi/180;
            x = x.*(ones(size(x))*cos(rot)) + ...
            y.*(ones(size(y))*sin(rot));

            %force inputs to be even
            if rem(length(x),2) == 1
                x = x(1:end-1,:);
                [m,n] = size(x);
            else
                [m,n] = size(x);
            end
            
            %spectra settings: shorter window/noverlap lengths = more
            %smoothing of resultant spectra
            x = detrend(x);
            nfft = 0.25*m;
            window = hanning(m,'periodic');
            noverlap = 0.7*length(window);
            
            bins = [1 5 10 15 20 25 35];
            Sx = NaN(m,n);
            F = NaN(m,n);
            for jj = bins
                [sx,f] = pwelch(x(:,jj),window,noverlap,nfft,fs); %compute PSD
                %locate data only between f = 1.2 - 10 Hz (Turbulent f range)
                fc = find(f >= 1.5 & f <= 15);
                sx = sx(fc);f = f(fc);
                ll = length(sx);
                Sx(1:ll,jj) = sx;F(1:ll,jj) = f;
            end
            clear f sx
            sx = Sx(~all(isnan(Sx),2),~all(isnan(Sx),1));
            f = F(~all(isnan(F),2),~all(isnan(F),1));
            clear Sx F
            
            %test plot
            [m,n] = size(sx);
            c = jet(n);
            sl = zeros(n,1);
            for jj = 1:n
                p = polyfit(log(f(:,jj)), log(sx(:,jj)),1);
                sl(jj) = p(1);
                yint(jj) = p(2);
                if ploton
                    figure(1)
                    y_hat=exp(p(1)*log(f(:,jj))+p(2));
                    k(1) = loglog(f(:,jj),y_hat,'--k'); hold on
                    g(jj) = loglog(f(:,jj),sx(:,jj),'+','Color',c(jj,:));
                end
            end
            if ploton
                legend([g(1),g(2),g(3),g(4),g(5),g(6),g(7)],'Bin 1','Bin 5','Bin 10','Bin 15','Bin 20','Bin 25','Bin 35')
                ylabel('S_z')
                xlabel('Frequency (Hz)')
                hold off
            end
            S.(fn{ii}).time(j,:) = VTA.(fn{ii}).time(ind(j)); %timestamp at the beginning of the two minute spectra
            S.(fn{ii}).rb(1,1:n) = VTA.(fn{ii}).rb(bins);
            S.(fn{ii}).hab(1,1:n) = hab - VTA.(fn{ii}).rb(bins);
            %Save spectral parameters for additional plots
            S.(fn{ii}).sp.f(j,:) = f(:,1);
            S.(fn{ii}).sp.bin1(j,:) = sx(:,1);
            S.(fn{ii}).sp.bin5(j,:) = sx(:,2);
            S.(fn{ii}).sp.bin10(j,:) = sx(:,3);
            S.(fn{ii}).sp.bin15(j,:) = sx(:,4);
            S.(fn{ii}).sp.bin20(j,:) = sx(:,5);
            S.(fn{ii}).sp.bin25(j,:) = sx(:,6);
            S.(fn{ii}).sp.bin35(j,:) = sx(:,7);
            S.(fn{ii}).yint(j,1:n) = yint; %intercept, p(2) from polyfit of spectra
            S.vpro1.hab(1,7) = 0; %for VP1, set the last rb to zero (because it's below the bed surface)
            S.(fn{ii}).sl(j,1:n) = sl;
        end
    end
    Spec.(fn{ii}) = S.(fn{ii}); %timestamp at the beginning of the two minute spectra
end
clearvars -except Spec ploton H datdir fname hn fn ifspec

%calculate polyfits (degree 3 or 4) for the slope data (for nice display)
fn = fieldnames(Spec);
for ii = 1:3                                %instruments
    [~,n] = size(Spec.(fn{ii}).sl);
    for j = 1:n                             %bin numbers
        sl = Spec.(fn{ii}).sl(:,j);
        nanid = find(isnan(sl));
        if length(nanid) > 0 %#ok<ISMT>
            %interp nans
            sl = cmgbridge(sl,40,200,200);
        end
        pw = 3;
        pf = polyfit((1:length(sl))',sl,pw);
        pv = polyval(pf,(1:length(sl))');
        count = 1;
        for jj = 1:5:length(sl)-1
            slstd(count) = std(sl(jj:jj+1));
            count = count+1;
        end
        slx = linspace(0,length(sl),length(slstd));
        sly = pv(1:5:end)';
        ebtime = Spec.(fn{ii}).time(1:5:end);
        %             pv(nanid) = NaN;                                %put nans back in wherever they were before calculating polyval
        Spec.(fn{ii}).slfit(:,j) = pv;
        Spec.(fn{ii}).err.slx(:,j) = slx;
        Spec.(fn{ii}).err.sly(:,j) = sly;
        Spec.(fn{ii}).err.slstd(:,j) = slstd;
        Spec.(fn{ii}).err.sltime(:,j) = ebtime;
        clear slstd ebtime sly slx pv
    end
end

% clearvars -except Spec n hn fn datdir ifspec
if ifspec
    prompt = 'Datafile all ready exists in the folder, overwrite? [y/n] ';
    result = input(prompt,'s');
    if strcmp(result,'y');
        save([datdir 'VTASlopeSpectra'],'Spec','-v7.3')
        disp('Datafile saved')
    elseif strcmp(result,'n');
    end
elseif ~ifspec
    save([datdir 'VTASlopeSpectra'],'Spec','-v7.3')
    disp('Datafile saved')
end
%% Plot routine
c = linspecer(n,'sequential');
f2 = figure(2);
set(f2,'PaperOrientation','portrait',...
    'position',[400 100   800 1000]);
set(gcf,'color','w','PaperPositionMode','auto')
%auto plotting with matrices! Now I feel extra smart...
pn = [3 2 1];
rc = 1;
ebspc = 1:2:14;
for ii = 1:3
    p(rc) = subplot(3,1,pn(rc));
    time = Spec.(fn{ii}).time;
    hab = Spec.(fn{ii}).hab;
    x = linspace(time(1),time(end),length(time));
    y = -5/3*ones(length(time),1);
    line(x,y,'Linewidth',3,'Color','k','LineStyle','--')
    hold on
    for j = 1:n
        r(j) = plot(time,Spec.(fn{ii}).slfit(:,j),...
            'Color',c(j,:),'Linewidth',1.5,...
            'DisplayName',sprintf('HAB = %.3fm',hab(j)));
        eb(j) = errorbar(Spec.(fn{ii}).err.sltime(ebspc(j):40:end,j),...
            Spec.(fn{ii}).err.sly(ebspc(j):40:end,j),...
            Spec.(fn{ii}).err.slstd(ebspc(j):40:end,j),'.',...
            'MarkerSize',1,'Color',c(j,:),'Linewidth',1.5);hold on
    end
    hold off
    set(gca,'Xlim',[time(1) time(end)],'YLim',[-3.5 1],'YTick',-3:1:1,...
        'LineWidth',1.5,'FontSize',16,'FontName','Cambria',...
        'TickDir','in','TickLength',[0.02 0.02],'box','on')
    
%     leg(rc) = legend(r);
    datetick('x','HH:MM','keepticks','keeplimits')
    %labeling
    if rc > 1
        set(gca,'XTickLabel',[])
    end
    if rc == 2
        ylabel('\bf\itPower Law Slope','FontName','Cambria')
        
    end
    if rc == 1
        xlabel('\bf\itTime on 14/03/2015','FontName','Cambria')
    end
    
    rc = rc+1;
end
% set(leg,'FontAngle','italic','FontWeight','Bold','LineWidth',1.5)
%block 1
set(p(1),'Position',[0.15 0.15 0.65 0.25])
set(p(2),'Position',[0.15 0.41 0.65 0.25])
set(p(3),'Position',[0.15 0.67 0.65 0.25])
% set(leg(1),'Position',[0.68 0.15 0.01 0.25])
% set(leg(2),'Position',[0.68 0.41 0.01 0.25])
% set(leg(3),'Position',[0.68 0.67 0.01 0.25])

figdir = 'c:\Users\bkn5\Projects\Mekong_W2015\Figures\Spectra&Waves\VerticalSpectra\';
export_fig([figdir 'VTASpecSlopesNoLeg'],'-pdf')
