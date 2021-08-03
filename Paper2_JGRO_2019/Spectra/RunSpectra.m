%Run spectral analysis: make a video of select frames.
%First, use default values for the spectra. We will play with these and
%figure out which values work best overall.
clear
savedatdir = 'C:\users\bkn5\Projects\Mekong_W2015\DataAnalysis\Spectra\Paper1\QCd\';
figdir = 'c:\Users\bkn5\Projects\Mekong_W2015\DataAnalysis\Spectra\Paper1\Spectra\SpecVPDepths\';
start = datenum(2015,03,07,14,24,00); %random sample
finish = datenum(2015,03,07,14,26,00);
name = ['VP2HTA_' datestr(start,'HHMMSS')];
inst = 'vpro2'; %'vpro2'; 'vpro3';
vph = 0.063; %height of VP off the bottom
colormap hot
load([savedatdir 'HTAday1Vels.mat'])
hab = vph-HTA.(inst).rb;
fs = HTA.(inst).sr;

%window the time series. For now, use a 10sec window and 1 second step
win = 10; %second
step = 0.1; %second
nsamp = win*fs;
intv = find(HTA.(inst).time >= start & HTA.(inst).time <= finish);
ind = intv(1):step*fs:intv(end);
ind = [ind intv(end)];

%create MPEG file; MPEG file settings
aviobj = VideoWriter([figdir name],'MPEG-4');
aviobj.FrameRate = 10;
aviobj.Quality = 100;
open(aviobj);
for i = 1:length(ind) %loop through steps
            if abs(length(HTA.(inst).z)-ind(i)) < nsamp %this deals with the remainder of the t-s when it is less than nsamp (only near the end of the t-s)
                idx = ind(i):length(HTA.(inst).z); 
                time = HTA.(inst).time(idx(1)); %get time stamp for plotting
                z = HTA.(inst).z(idx,:);
            else
                idx = ind(i):ind(i)+nsamp; %steps include number of samples equal to the window size * the sampling freq
                time = HTA.(inst).time(idx(1)); %get time stamp for plotting
                z = HTA.(inst).z(idx,:);
            end
    
    if isempty(z)
        continue
    else
        if rem(length(z),2) == 1
            z = z(1:end-1,:);
            [m,n] = size(z); %force inputs to be even
        else
            [m,n] = size(z);
        end
        z = detrend(z);
        nfft = 0.25*m;
        window = hanning(m,'periodic');
        noverlap = 0.7*length(window);
        minf = 0;
        maxf = 10;
        
        Sz = NaN(m,n);
        F = NaN(m,n);
        for k = 1:n
            [sz,f] = pwelch(z(:,k),window,noverlap,nfft,fs); %compute PSD
            fr = find(f >= minf & f <= maxf); %set frequency cutoff
            sz = sz(fr);f = f(fr);
            ll = length(sz);
            Sz(1:ll,k) = sz;F(1:ll,k) = f;
        end
        clear f sz
        sz = Sz(~all(isnan(Sz),2),~all(isnan(Sz),1));
        f = F(~all(isnan(F),2),~all(isnan(F),1));
        
        id = find(hab<0);
        if ~isempty(id)
            id = id(end);
        end
        f1 = figure(1);
        set(f1,'PaperOrientation','portrait',...
            'position',[400 100   800 600]);
        set(gcf,'color','w','PaperPositionMode','auto')
        imagesc(f(:,1),hab,flipud(rot90(log10(sz))));
        caxis([-4.5 -3])
        set(gca,'Xlim',[0 10],'XTick',0:1:10,'Ylim',[min(hab) max(hab)],'YDir','normal')
        line(linspace(0,10,length(f(:,1))),zeros(length(f(:,1)),1),'LineWidth',2,'Color','w')
        xlabel('Frequency (Hz)')
        ylabel('Measurement height above bottom (m)')
        title(['Spectra window length ' num2str(win) 's at time ' datestr(time,'HH:MM:SS:FFF')])
        cb = colorbar;
        ylabel(cb,'Log_1_0(S_z)')
        %write frame to avi file
        frame = getframe(f1);
        writeVideo(aviobj,frame);
        %             drawnow
        %             frame = getframe(1);
        %             im = frame2im(frame);
        %             [imind,cm] = rgb2ind(im,256);
        %             if i == 1;
        %                 imwrite(imind,cm,[figdir name],'gif', 'Loopcount',inf,'DelayTime',0.033);
        %             else
        %                 imwrite(imind,cm,[figdir name],'gif','WriteMode','append','DelayTime',0.033);
        %             end
    end
end
close(aviobj)
close all