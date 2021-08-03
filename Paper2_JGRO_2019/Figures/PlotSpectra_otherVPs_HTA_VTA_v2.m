%Plot cross-shore velocity spectra from the other instruments other than
%the St # peak instreuments.
clear

figdir = 'f:\GradSchool\DataAnalysis\Paper2\WorkingFigures\Spectra\';
datdir = 'd:\Projects\Mekong_W2015\DataAnalysis\Paper3\VPs\';
load('d:\Projects\Mekong_W2015\DataAnalysis\Paper2\Spectra\StSpect_rawData.mat')
fn = fieldnames(data);
txt = {'x = -10 cm','x = 10 cm';'x = -10 cm','x = 20 cm';'z/h_c = 0.6','z/h_c = 1.25'};
%%%Plot Routine

f1 = figure(1);
set(f1,'PaperOrientation','portrait',...
    'position',[400 100   1000   400],...
    'renderer','painters');
line = {'-';'--';'-.'};
heading = [340 240 240]; %headings are in order! (HTA1 - 3, VTA)
vps = [1 2; 1 3; 2 3];
sp = zeros(3,1);
for i = 1:3
    dn = fieldnames(data.(fn{i}));
    if i == 1
        bins = 1:5;
    else
        bins = 9:23;
    end
    %Spectra settings
    fs = 50;
    nwin = fs*600;
    swin = fs*50; %20 second averaging window (to smooth)
    DOF = round((nwin/swin)*2);
    if i == 1
        disp(['50% Hamming windowed spectra with ' sprintf('%0.0f',DOF) ' degrees of freedom'])
    end
    clr = [0 0 0;0.5 0.5 0.5];
    b = zeros(2,1);
    lns = {'-';'-'};
    for ii = 1:2
        j = vps(i,ii);
        sp(i) = subplot(1,3,i);
        x = data.(fn{i}).(dn{j}).x;
        y = data.(fn{i}).(dn{j}).y;
%         rot = (pi*heading(i))/180;
%         T = [cos(rot) -sin(rot);...
%             sin(rot) cos(rot)];
%         vels = [x y];
%         V = vels*T';
%         x = V(:,1);y = V(:,2);
        
        %calculate psd
        [Cuu,F,cf] = pwelch(detrend(x),hanning(swin),swin*0.5,nwin,fs,'confidencelevel',0.95);
        Cuu(F < 0.05) = NaN;
        %confidence interval
        L = nanmean(10*log10(Cuu)-10*log10(cf(:,1)));
        U = nanmean(10*log10(cf(:,2))-10*log10(Cuu));
        
        %plot 5/3 slope on figures
        if i == 3
            int = 0.0015;
            text(4,2E-4,'^-^5^/^3'),hold on
        else
            int = 0.005;
            text(4,6E-4,'^-^5^/^3'),hold on
        end
        hold on
        if ii == 1 && i ~=4
            xs = linspace(2,45,length(Cuu));
            ys = int.*(xs.^(-5/3));
            plot(xs,ys,'--k','LineWidth',1.5);
        end
        b(ii) = plot(F,Cuu,...
            'linestyle',lns{ii},...
            'Color',clr(ii,:),...
            'LineWidth',1.5);
        if ii == 1
            errorbar(4,1E-2,L*1E-2,U*1E-2,...
                'marker','.',...
                'linewidth',1.5,...
                'color','k')
            text(5,1E-2,'\bfx 10')
        end
        set(gca,'yscale','log','xscale','log')
    end
    leg = legend(b,txt{i,:});
    set(leg,...
        'box','off',...
        'location','southwest')
    
    clear dat
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
ylabel(sp(1),'S_{uu} (m^2/s^2/Hz)')
title(sp(1),'HTA, z/h_c = 0.03')
title(sp(2),'HTA, z/h_c = 0.33')
title(sp(3),'VTA')
% export_fig([figdir 'HTA_VTAotherVPs_v2'],'-pdf')
