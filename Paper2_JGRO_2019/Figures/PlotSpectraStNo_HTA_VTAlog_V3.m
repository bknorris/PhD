%Plot cross-shore velocity spectra from the 3 VPs of the HTA and the 3 VPs
%of the VTA experiment. Input time series are 30 minutes of data during the
%individual experiments.
clear, close all

figdir = 'f:\GradSchool\DataAnalysis\Paper2\WorkingFigures\Spectra\';
load('d:\Projects\Mekong_W2015\DataAnalysis\Paper2\Spectra\StSpect_rawData.mat')    
fn = fieldnames(data);
txt = {'x = 20 cm';'x = 10 cm';''};
%%%Plot Routine

f1 = figure(1);
set(f1,'PaperOrientation','portrait',...
    'position',[400 100   1000   400],...
    'renderer','painters');
line = {'-';'--';'-.'};
heading = [340 240 240]; %headings are in order! (HTA1 - 3, VTA)
vps = [3 2 1];
sp = zeros(3,1);
for i = 1:3
    for j = 1:3
        if j == vps(i)
            dn = fieldnames(data.(fn{i}));
            %Spectra settings
            fs = 50;
            nwin = fs*600;
            swin = fs*110; %20 second averaging window (to smooth)
            DOF = round((nwin/swin)*2);
            disp(['50% Hamming windowed spectra with ' sprintf('%0.0f',DOF) ' degrees of freedom'])
            
            sp(i) = subplot(1,3,i);
            x = data.(fn{i}).(dn{vps(i)}).x;
            y = data.(fn{i}).(dn{vps(i)}).y;

            %calculate psd
            [Cuu,F,cfu] = pwelch(detrend(x),hanning(swin),swin*0.5,nwin,fs,'confidencelevel',0.95);
            Cuu(F < 0.05) = NaN;
            %confidence interval
            L = nanmean(10*log10(Cuu)-10*log10(cfu(:,1)));
            U = nanmean(10*log10(cfu(:,2))-10*log10(Cuu));

            %plot 5/3 slope on figures
            if i == 3
                int = 0.0015;
                text(4,2E-4,'^-^5^/^3'),hold on
            else
                int = 0.005;
                text(4,6E-4,'^-^5^/^3'),hold on
            end
            xs = linspace(2,45,length(Cuu));
            ys = int.*(xs.^(-5/3));
            plot(xs,ys,'-.k','LineWidth',1.5);hold on
            loglog(F,Cuu,...
                'Color','k',...
                'LineWidth',2)
            errorbar(4,1E-2,L*1E-2,U*1E-2,...
                'marker','.',...
                'linewidth',1.5,...
                'color','k')
            text(5,1E-2,'\bfx 10')
            text(0.15,2E-5,txt{i},'FontSize',14)
            set(gca,'yscale','log','xscale','log')
        end
    end
end
%plot adjustments
set(sp,'xlim',[0.1 10],'ylim',[1E-5 1E-1])
set(sp(1),'position',[0.1 0.15 0.26 0.7])
set(sp(2),'position',[0.405 0.15 0.26 0.7],...
    'yticklabel',[])
set(sp(3),'position',[0.71 0.15 0.26 0.7],...
    'yticklabel',[])
xlabel(sp(2),'f (Hz)')
ylabel(sp(1),'S_{uu} (m^2/s^2/Hz)')
title(sp(1),'HTA, z/h_c = 0.03')
title(sp(2),'HTA, z/h_c = 0.33')
title(sp(3),'VTA, z/h_c = 0.04')
prettyfigures('text',13,'labels',14,'box',1)
% export_fig([figdir 'HTA_VTApsd_lsV2_1'],'-pdf')
        
