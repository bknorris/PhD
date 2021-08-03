%Create evolutionary Spectra plots for all 3 days, all 3 instruments. Use a
%variety of bins to look for something interesting.

datdir = 'c:\Users\bkn5\Projects\Mekong_W2015\DataAnalysis\Spectra\Paper1\QCd\';
figdir = 'c:\Users\bkn5\Projects\Mekong_W2015\Figures\Spectra&Waves\VerticalSpectra\';
load([datdir 'HTASlopeSpectra.mat'])
hn = {'day1';'day2';'day3'};
fn = fieldnames(Spec.day1);
sfn = fieldnames(Spec.day1.vpro1.sp);
plotn = [1 2 3;4 5 6;7 8 9];
for i = 1:3
    for ii = 1:3
        f(plotn(i,ii)) = figure(plotn(i,ii));
        set(f(plotn(i,ii)),'PaperOrientation','portrait',...
            'position',[400 100   800 1000]);
        set(gcf,'color','w','PaperPositionMode','auto')
        colormap hot
        for j = 1:length(sfn)
            if strcmp(sfn{j},'f')
                continue
            else
                s(j) = subplot(length(sfn),1,j);
                C = log10(Spec.(hn{i}).(fn{ii}).sp.(sfn{j}));
                y = Spec.(hn{i}).(fn{ii}).sp.f(1,:);
                x = Spec.(hn{i}).(fn{ii}).time;
                p(j) = imagesc(x,y',C');
                if i == 1; 
                    tname = {'dummy';'bin1';'bin2';'bin3';'bin4';'bin5';'bin6';'bin7'};
                else
                    tname = sfn;
                end
                title(tname{j},'FontSize',14)
                set(gca,'XLim',[x(1) x(end)],'YLim',[1.5 10],'YTick',2:2:10,'FontSize',11,'LineWidth',1.5)
                if j == 5
                    ylabel('\bf\itFrequency (Hz)','FontSize',14)
                end
                datetick('x','HH:MM','keepticks','keeplimits')
                if j < 8
                    set(gca,'XTickLabel',[])
                end
                if j == 8
                    xlabel(['\bf\itTime on ' datestr(x(1),'dd/mm/yyyy')],'FontSize',14)
                end
                caxis([-6 -3])
            end
        end
        cb = colorbar;
        set(s(2),'position',[0.08 0.85 0.8 0.1])
        set(s(3),'position',[0.08 0.72 0.8 0.1])
        set(s(4),'position',[0.08 0.59 0.8 0.1])
        set(s(5),'position',[0.08 0.46 0.8 0.1])
        set(s(6),'position',[0.08 0.33 0.8 0.1])
        set(s(7),'position',[0.08 0.20 0.8 0.1])
        set(s(8),'position',[0.08 0.07 0.8 0.1])
        ylabel(cb,'\bf\itLog_1_0(m^-^2s^-^1)','FontSize',14)
        set(cb,'position',[0.9 0.07 0.02 0.88],'LineWidth',1.5,'FontSize',11)
        
        %save figure
        name = ['VertSpEv' hn{i} fn{ii}];
        export_fig([figdir name],'-png')
    end
end
close all
datdir = 'c:\Users\bkn5\Projects\Mekong_W2015\DataAnalysis\Spectra\Paper1\';
load([datdir 'VTASlopeSpectra.mat'])
fn = fieldnames(Spec);
sfn = fieldnames(Spec.vpro1.sp);
for i = 1:3
    for j = 1:length(sfn)
        f(i) = figure(i);
        set(f(i),'PaperOrientation','portrait',...
            'position',[400 100   800 1000]);
        set(gcf,'color','w','PaperPositionMode','auto')
        colormap hot
        if strcmp(sfn{j},'f')
            continue
        else
            s(j) = subplot(length(sfn),1,j);
            C = log10(Spec.(fn{i}).sp.(sfn{j}));
            y = Spec.(fn{i}).sp.f(1,:);
            x = Spec.(fn{i}).time;
            p(j) = imagesc(x,y',C');
            title(sfn{j},'FontSize',14)
            set(gca,'XLim',[x(1) x(end)],'YLim',[1.5 10],'YTick',2:2:10,'FontSize',11,'LineWidth',1.5)
            if j == 5
                ylabel('\bf\itFrequency (Hz)','FontSize',14)
            end
            datetick('x','HH:MM','keepticks','keeplimits')
            if j < 8
                set(gca,'XTickLabel',[])
            end
            if j == 8
                xlabel(['\bf\itTime on ' datestr(x(1),'dd/mm/yyyy')],'FontSize',14)
            end
            caxis([-6 -4])
        end
    end
    cb = colorbar;
    set(s(2),'position',[0.08 0.85 0.8 0.1])
    set(s(3),'position',[0.08 0.72 0.8 0.1])
    set(s(4),'position',[0.08 0.59 0.8 0.1])
    set(s(5),'position',[0.08 0.46 0.8 0.1])
    set(s(6),'position',[0.08 0.33 0.8 0.1])
    set(s(7),'position',[0.08 0.20 0.8 0.1])
    set(s(8),'position',[0.08 0.07 0.8 0.1])
    ylabel(cb,'\bf\itLog_1_0(m^-^2s^-^1)','FontSize',14)
    set(cb,'position',[0.9 0.07 0.02 0.88],'LineWidth',1.5,'FontSize',11)
    
    %save figure
    name = ['VertSpEv' 'VTA' fn{i}];
    export_fig([figdir name],'-png')
end
close all
    
