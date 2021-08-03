%Calculate the slope of velocity spectra from VecPros during the HTA and
%VTA. Calculate spectra at set heights: Bins 1,5,10,15,20,25,35 for each
%VecPro during the experiments. Capture the slope of the spectra from 1.5Hz
%to 10Hz, and save it to a variable.
clear
close all
ploton = 0;
datdir = 'c:\Users\bkn5\Projects\Mekong_W2015\DataAnalysis\Paper1\';
fdir = dir([datdir 'HTAday*Vels.mat']);
fname = {fdir.name};
Spec = struct();
%define instrument heights from deployment logs
H.day1.vpro1 = 0.062;H.day1.vpro2 = 0.063;H.day1.vpro3 = 0.063;
H.day2.vpro1 = 0.242;H.day2.vpro2 = 0.271;H.day2.vpro3 = 0.240;
H.day3.vpro1 = 0.550;H.day3.vpro2 = 0.555;H.day3.vpro3 = 0.535;
%Check for the existence of the datafile to conditionally run the program
ifspec = exist([datdir 'HTASlopeSpectra.mat'],'file');
for i = 1:3
    disp(['Loading ' fname{i}])
    load([datdir fname{i}])
    hn = fieldnames(H);
    S = struct();
    %loop through the instruments
    for ii = 2:4
        fn = fieldnames(HTA);
        disp(['Running spectral analysis on ' fn{ii}])
        hab = H.(hn{i}).(fn{ii});
        fs = HTA.(fn{ii}).sr;
        
        %window the time series. For now, use a 30sec window and 10 second step
        win = 30; %second
        step = 10; %second
        nsamp = win*fs;
        ind = 1:step*fs:length(HTA.(fn{ii}).z1);
        for j = 1:length(ind)
            if abs(length(HTA.(fn{ii}).z1)-ind(j)) < nsamp %skip the last few indexes approaching the end of the t-s
                continue
            else
                idx = ind(j):ind(j)+nsamp; %steps include number of samples equal to the window size * the sampling freq
                x = HTA.(fn{ii}).x(idx,:);
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
                minf = 3.7;
                maxf = 15; %determined using the frozen turb theory/analysis of spectral slopes
                
                if i == 1 %i.e. on the first day, use bins 1-7
                    bins = 1:7;
                else
                    bins = [1 5 10 15 20 25 35];
                end
                Sz = NaN(m,n);
                F = NaN(m,n);
                for jj = bins
                    [sz,f] = pwelch(x(:,jj),window,noverlap,nfft,fs); %compute PSD
                    %locate data only between f = 6 - 15 Hz (Turbulent f range?)
                    fc = find(f >= minf & f <= maxf);
                    sz = sz(fc);f = f(fc);
                    ll = length(sz);
                    Sz(1:ll,jj) = sz;F(1:ll,jj) = f;
                end
                clear f sz
                sz = Sz(~all(isnan(Sz),2),~all(isnan(Sz),1));
                f = F(~all(isnan(F),2),~all(isnan(F),1));
                clear Sz F
                
                %test plot
                [m,n] = size(sz);
                c = jet(n);
                sl = zeros(n,1);
                yint = zeros(n,1);
                for jj = 1:n
                    p = polyfit(log(f(:,jj)), log(sz(:,jj)),1);
                    sl(jj) = p(1);
                    yint(jj) = p(2);
                    if ploton
                        figure(1)
                        y_hat=exp(p(1)*log(f(:,jj))+p(2));
                        k(1) = loglog(f(:,jj),y_hat,'--k'); hold on
                        g(jj) = loglog(f(:,jj),sz(:,jj),'+','Color',c(jj,:));
                    end
                end
                if ploton
                    legend([g(1),g(2),g(3),g(4),g(5),g(6),g(7)],'Bin 1','Bin 5','Bin 10','Bin 15','Bin 20','Bin 25','Bin 35')
                    ylabel('S_z')
                    xlabel('Frequency (Hz)')
                    hold off
                end
                S.(fn{ii}).time(j,:) = HTA.(fn{ii}).time(ind(j)); %timestamp at the beginning of the two minute spectra
                S.(fn{ii}).rb(1,1:n) = HTA.(fn{ii}).rb(bins);
                S.(fn{ii}).hab(1,1:n) = hab - HTA.(fn{ii}).rb(bins);
                %Save spectral parameters for additional plots
                S.(fn{ii}).sp.f(j,:) = f(:,1);
                S.(fn{ii}).sp.bin1(j,:) = sz(:,1);
                S.(fn{ii}).sp.bin5(j,:) = sz(:,2);
                S.(fn{ii}).sp.bin10(j,:) = sz(:,3);
                S.(fn{ii}).sp.bin15(j,:) = sz(:,4);
                S.(fn{ii}).sp.bin20(j,:) = sz(:,5);
                S.(fn{ii}).sp.bin25(j,:) = sz(:,6);
                S.(fn{ii}).sp.bin35(j,:) = sz(:,7);
%                 if i == 1
%                     S.(fn{ii}).hab(1,7) = 0; %on day1, set the last rb to zero (because it's below the bed surface)
%                 end
                S.(fn{ii}).sl(j,1:n) = sl;
                S.(fn{ii}).yint(j,1:n) = yint; %intercept, p(2) from polyfit of spectra
            end
        end
        Spec.(hn{i}).(fn{ii}) = S.(fn{ii}); %timestamp at the beginning of the two minute spectra
    end
    clearvars -except Spec ploton H datdir fname hn fn ifspec
end

hn = fieldnames(Spec);fn = fieldnames(Spec.day1);
%calculate polyfits (degree 3 or 4) for the slope data (for nice display)
for i = 1:3                                     %days
    for ii = 1:3                                %instruments
        [~,n] = size(Spec.(hn{i}).(fn{ii}).sl);
        for j = 1:n                             %bin numbers
            sl = Spec.(hn{i}).(fn{ii}).sl(:,j);
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
            ebtime = Spec.(hn{i}).(fn{ii}).time(1:5:end);
%             pv(nanid) = NaN;                                %put nans back in wherever they were before calculating polyval
            Spec.(hn{i}).(fn{ii}).slfit(:,j) = pv;
            Spec.(hn{i}).(fn{ii}).err.slx(:,j) = slx;
            Spec.(hn{i}).(fn{ii}).err.sly(:,j) = sly;
            Spec.(hn{i}).(fn{ii}).err.slstd(:,j) = slstd;
            Spec.(hn{i}).(fn{ii}).err.sltime(:,j) = ebtime;
            clear slstd ebtime sly slx pv
        end
    end
end
clearvars -except Spec n hn fn datdir ifspec
if ifspec
    prompt = 'Datafile all ready exists in the folder, overwrite? [y/n] ';
    result = input(prompt,'s');
    if strcmp(result,'y');
        save([datdir 'HTASlopeSpectra'],'Spec','-v7.3')
        disp('Datafile saved')
    elseif strcmp(result,'n');
    end
elseif ~ifspec
    save([datdir 'HTASlopeSpectra'],'Spec','-v7.3')
    disp('Datafile saved')
end

%% Plot routine
c = linspecer(n,'sequential');
f2 = figure(2);
set(f2,'PaperOrientation','portrait',...
    'position',[400 100   1400 1000]);
set(gcf,'color','w','PaperPositionMode','auto')
%auto plotting with matrices! Now I feel extra smart...
pn = [7 8 9;4 5 6;1 2 3]; %subplot number
rc = 1;
cc = 1;
ebspc = 1:2:14;
for i = 1:3
    for ii = 1:3
        p(pn(rc,cc)) = subplot(3,3,pn(rc,cc));
        time = Spec.(hn{i}).(fn{ii}).time;
        hab = Spec.(hn{i}).(fn{ii}).hab;
        x = linspace(time(1),time(end),length(time));
        y = -5/3*ones(length(time),1);
        line(x,y,'Linewidth',3,'Color','k','LineStyle','--')
        hold on
        for j = 1:n
            r(j) = plot(time,Spec.(hn{i}).(fn{ii}).slfit(:,j),...
                'Color',c(j,:),'Linewidth',1.5,...
                'DisplayName',sprintf('HAB = %.3fm',hab(j)));
            eb(j) = errorbar(Spec.(hn{i}).(fn{ii}).err.sltime(ebspc(j):40:end,j),...
                Spec.(hn{i}).(fn{ii}).err.sly(ebspc(j):40:end,j),...
                Spec.(hn{i}).(fn{ii}).err.slstd(ebspc(j):40:end,j),'.',...
                'MarkerSize',1,'Color',c(j,:),'Linewidth',1.5);hold on
        end
        hold off
        set(gca,'Xlim',[time(1) time(end)],'YLim',[-3 0],'YTick',-3:1:0,...
            'LineWidth',1.5,'FontSize',16,'FontName','Cambria',...
            'TickDir','in','TickLength',[0.02 0.02],'box','on')
        %legend only for third block in the row
        if cc == 3
%             leg(i) = legend(r);
        end
        %labeling
        if ii >= 2
            set(gca,'YTickLabel',[])
        end
        if cc == 1 && rc == 2
            ylabel('\bf\itPower Law Slope','FontName','Cambria')
        end
        if cc == 2 && rc == 3
            xlabel('\bf\itTime on 10/03/2015','FontName','Cambria')
        elseif cc == 2 && rc == 2
            xlabel('\bf\itTime on 08/03/2015','FontName','Cambria')
        elseif cc == 2 && rc == 1
            xlabel('\bf\itTime on 07/03/2015','FontName','Cambria')
        end
        datetick('x','HH:MM','keepticks','keeplimits')
        cc = cc+1;
    end
    rc = rc+1;
    cc = 1;
end
% set(leg,'FontAngle','italic','FontWeight','Bold','LineWidth',1.5)
%block 1
set(p(7),'Position',[0.08 0.07 0.22 0.25])
set(p(8),'Position',[0.32 0.07 0.22 0.25])
set(p(9),'Position',[0.56 0.07 0.22 0.25])
% set(leg(1),'Position',[0.88 0.07 0.01 0.25])
%block 2
set(p(4),'Position',[0.08 0.4 0.22 0.25])
set(p(5),'Position',[0.32 0.4 0.22 0.25])
set(p(6),'Position',[0.56 0.4 0.22 0.25])
% set(leg(2),'Position',[0.88 0.4 0.01 0.25])
%block 3
set(p(1),'Position',[0.08 0.72 0.22 0.25])
set(p(2),'Position',[0.32 0.72 0.22 0.25])
set(p(3),'Position',[0.56 0.72 0.22 0.25])
% set(leg(3),'Position',[0.88 0.72 0.01 0.25])

figdir = 'c:\Users\bkn5\Projects\Mekong_W2015\Figures\Spectra&Waves\VerticalSpectra\';
export_fig([figdir 'HTASpecSlopesNoLeg'],'-pdf')


