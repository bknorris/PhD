%Plot cross-shore velocity spectra from the other instruments other than 
%the St # peak instreuments. 
clear

figdir = 'd:\Projects\Mekong_W2015\Figures\Paper2\';
datdir = 'd:\Projects\Mekong_W2015\DataAnalysis\Paper1\';
files = {'HTA_1Vels.mat';'HTA_2Vels.mat';'VTA_2vp2Vels.mat';'VTA_2vp3Vels.mat'};
txt = {'x = -10 cm','x = 10 cm';'x = -10 cm','x = 20 cm';'z/h_c = 0.6','z/h_c = 1.25'};
%%%Plot Routine

f1 = figure(1);
set(f1,'PaperOrientation','portrait',...
    'position',[400 100   1000   400],...
    'renderer','painters');
line = {'-';'--';'-.'};
heading = [20 20 96 96];
vps = [1 2; 1 3; 1 0;0 1];
sp = zeros(3,1);
cc = 1;
for i = 1:4
    disp(['Loading ' files{i}])
    load([datdir files{i}])
    fn = fieldnames(dat);
    if i == 4
        cc = 3;
    end
    if i > 1
        bins = 1:5;
    else
        bins = 9:23;
    end
    %Spectra settings
    fs = 50;
    win = 300; %seconds (5 minutes)
    nwin = fs*win;
    swin = fs*20; %20 second averaging window (to smooth)
    DOF = round((nwin/swin)*2);
    if i == 1
        disp(['50% Hamming windowed spectra with ' sprintf('%0.0f',DOF) ' degrees of freedom'])
    end
    clr = [0 0 0;0.5 0.5 0.5];
    b = zeros(2,1);
    lns = {'--';'-'};
    for ii = 1:2
        if vps(i,ii) == 0
            continue
        else
            j = vps(i,ii);
        end
        sp(cc) = subplot(1,3,cc);
        x = mean(dat.(fn{j}).x(:,bins),2);
        y = mean(dat.(fn{j}).y(:,bins),2);
        rot = (pi*heading(i))/180;
        T = [cos(rot) -sin(rot);...
            sin(rot) cos(rot)];
        vels = [x y];
        V = vels*T';
        x = V(:,1);y = V(:,2);
        
        %calculate psd
        x = x(1:nwin*6);
        [Cuu,F,cf] = pwelch(detrend(x),hanning(swin),swin*0.7,nwin,fs,'confidencelevel',0.95);
        Cuu(F < 0.05) = NaN;
        cf(F < 0.05,:) = [];
        ff = F;
        ff(F < 0.05,:) = [];
        %plot 5/3 slope on figures
        if i > 2
            int = 0.0015;
            text(3,5E-4,'^-^5^/^3'),hold on
        else
            int = 0.005;
            text(3,1E-3,'^-^5^/^3'),hold on
        end
        hold on
        if ii == 1 && i ~=4
            xs = linspace(2,45,length(Cuu));
            ys = int.*(xs.^(-5/3));
            plot(xs,ys,'Color','r','LineWidth',1.5);
        end
        xx = [ff;flipud(ff)];
        yy = [cf(:,1);flipud(cf(:,2))];
        fill(xx,yy,[0.8 0.8 0.8],...
            'EdgeColor','none'), hold on
        b(ii) = plot(F,Cuu,...
            'linestyle',lns{ii},...
            'Color',clr(ii,:),...
            'LineWidth',1.5);
        set(gca,'yscale','log','xscale','log')    
    end
    try
        leg = legend(b,txt{cc,:});set(leg,'box','off')
    catch
        continue
    end
    cc = cc+1;
    clear dat b
end

%plot adjustments
prettyfigures('text',13,'labels',14,'box',1)
set(sp,'xlim',[0.1 10])
set(sp,'ylim',[1E-6 1E-1])
set(sp(1),'position',[0.1 0.15 0.26 0.7])
set(sp(2),'position',[0.41 0.15 0.26 0.7],...
    'yticklabel',[])
set(sp(3),'position',[0.72 0.15 0.26 0.7],...
    'yticklabel',[])
xlabel(sp(2),'f (Hz)')
ylabel(sp(1),'Spectral Density (m^2s^-^1)')
title(sp(1),'HTA, z/h_c = 0.03')
title(sp(2),'HTA, z/h_c = 0.33')
title(sp(3),'VTA')
export_fig([figdir 'HTA_VTAotherVPs'],'-pdf')
