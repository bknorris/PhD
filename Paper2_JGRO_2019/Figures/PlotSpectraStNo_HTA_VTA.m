%Plot cross-shore velocity spectra from the 3 VPs of the HTA and the 3 VPs
%of the VTA experiment. Input time series are 30 minutes of data during the
%individual experiments.
clear

ddir = 'd:\Projects\Mekong_W2015\DataAnalysis\Paper1\Eps_Vel_Spectra\';
files = dir([ddir '*30min.mat']);
files = {files.name};
days = {'day1';'day2';'day3';'day4'};
%%%Plot Routine

f1 = figure(1);
set(f1,'PaperOrientation','portrait',...
    'position',[400 100   1000   700]);
cl = [0.1 0.1 0.1;0.5 0.5 0.5;0.7 0.7 0.7];
line = {'-';'--';'-.'};
pl = [5 3 1]; %HTA
pl2 = [6 4 2]; %VTA
heading = [20 20 20 96];
sp = zeros(6,1);
data = struct();
for i = 1:length(files)
    load([ddir files{i}])
    fn = fieldnames(dat);
    if i == 1
        bins = 1:5;
    else
        bins = 9:23;
    end
    %Spectra settings
    fs = 50;
    win = 900; %seconds (5 minutes)
    nwin = fs*win;
    swin = fs*60; %30 second averaging window (to smooth)
    DOF = round((nwin/swin)*2);
    if i == 1
        disp(['50% Hamming windowed spectra with ' sprintf('%0.0f',DOF) ' degrees of freedom'])
    end
    for ii = 1:3
        if i < 4
            sp(i) = subplot(3,2,pl(i));
            cc = flipud(cl);
        else
            sp(i+3) = subplot(3,2,pl2(ii));
            cc = cl;
        end
        x = mean(dat.(fn{ii}).x(:,bins),2);
        y = mean(dat.(fn{ii}).y(:,bins),2);
        rot = (pi*heading(i))/180;
        T = [cos(rot) -sin(rot);...
            sin(rot) cos(rot)];
        vels = [x y];
        V = vels*T';
        x = V(:,1);y = V(:,2);
        %calculate RMS velocities (Luhar et al. 2013)
        Ec = (1/length(y))*sum(y);Nc = (1/length(x))*sum(x);
        Ewrms = sqrt((1/length(y))*sum((y-Ec).^2));
        Nwrms = sqrt((1/length(x))*sum((x-Nc).^2));
        Uc = sqrt(Ec^2+Nc^2);Uwrms = sqrt(Ewrms^2+Nwrms^2);
        Uw = sqrt(2)*Uwrms;
        %calculate Eulerian time scale (autocorr of x)
        [acor,~] = xcorr(x,'coeff');
        acor = acor(length(x):end);
        id = find(acor <= 0,1,'first');
        tu = trapz(acor(1:id-1))/fs;
        %save data to structure
        data.(days{i}).(fn{ii}).Uc = Uc;
        data.(days{i}).(fn{ii}).Uw = Uw;
        data.(days{i}).(fn{ii}).tu = tu;
        %calculate psd
        [Cuu,F] = pwelch(x,hanning(nwin),swin*0.5,swin,fs);
        Cuu = my_running_median(Cuu,5);
        Cuu = smooth(Cuu,12);
        Cuu(F < 0.1) = NaN;
        loglog(F,Cuu,...
            'Color',cc(ii,:),...
            'LineWidth',1.5), hold on
    end
end
%plot 5/3 slope on figures
for i = 1:6
    subplot(3,2,i)
    int = 0.0022;
    xs = linspace(2,45,length(Cuu));
    ys = int.*(xs.^(-5/3));
    loglog(xs,ys,'Color','r','LineWidth',1.5);
    text(3.5,0.0005,'\itf ^-^5^/^3','FontName','Arial',...
        'FontSize',14)
end

    

        
        
