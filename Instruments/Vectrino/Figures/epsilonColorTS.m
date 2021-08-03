%make color plots of TKE dissipation rate for a few select deployments, 
%F2F, F2F2, F2F3, FSS2, and DPS
clear
files = {'F2F_1TKE.mat';'FSS2TKE.mat';'DPSTKE.mat'};
z = [108 208 308; 150 82 150; 150 82 150;84 179 279];
datdir1 = 'd:\Projects\Mekong_F2014\DataAnalysis\Paper2\TKE\';
fname = dir(datdir1);fname = {fname.name};

savefigdir = 'd:\Projects\Mekong_F2014\Figures\Paper2\';
for i = 1:length(fname)
    isfile = strcmp(fname{i},files);
    if any(isfile) == 1
        name = regexprep(fname{i},'TKE.mat','');
        load([datdir1 fname{i}])
        
        f1 = figure(i);
        set(f1,'PaperOrientation','portrait',...
            'position',[400 200   1200   500]);
        set(gcf,'color','w','PaperPositionMode','auto')
        c = 1;
        fn = fieldnames(Stat);
        for ii = 1:3
            E = (Stat.(fn{ii}).beam1.E+Stat.(fn{ii}).beam3.E)./2;
            time = Stat.(fn{ii}).time;
            rb =  [0.0400;0.0410;0.0419;0.0429;0.0438;0.0448;0.0457;0.0467;0.0477;0.0486;0.0496;0.0505;...
                0.0515;0.0524;0.0534;0.0543;0.0553;0.0563;0.0572;0.0582;0.0591;0.0601;0.0610;0.0620;...
                0.0630;0.0639;0.0649;0.0658;0.0668;0.0677;0.0687;0.0697;0.0706;0.0716;0.0725];            
            if strcmp(name,'DPS') || strcmp(name,'FSS2')
                h = rb;
            else
                h = (z(c,ii)./1000)-rb;
            end
            %plot
            p(ii) = subplot(1,3,ii);
            imagesc(time,h,E)
            [~,~,~,hr,mi,~] = datevec(time(end)-time(1));
            tint = (hr*60) + mi;                                                    %calculate minutes elapsed between start/stop time of experiment
            tstep = datenum(0,0,0,0,floor(tint/4),0);
            set(gca,'Xtick',time(1):tstep:time(end))
            datetick('x','HH:MM','keepticks','keeplimits')
            if ii == 2
                xlabel(['Time on ' datestr(time(1),'dd/mm/yyyy')],'FontSize',14,'FontName','Cambria')
                set(gca,'Yticklabel',[])
            end
            if ii == 1 && strcmp(name,'DPS') || ii == 1 && strcmp(name,'FSS2')
                ylabel('Distance from sensor (m)','FontSize',14,'FontName','Cambria')
            elseif ii == 1 && strcmp(name,'F2F_1')
                ylabel('Height above bed (m)','FontSize',14,'FontName','Cambria')
            end
            if ii == 3
                set(gca,'Yticklabel',[])
            end
            if c < length(z)
                c = c+1;
            end
        end
        %global plot adjustments
        set(p(1),'position',[0.08 0.12 0.25 0.76])
        set(p(2),'position',[0.36 0.12 0.25 0.76])
        set(p(3),'position',[0.64 0.12 0.25 0.76])
        cb = colorbar;
        ylabel(cb,'Dissipation Rate (m^2s^-^3)','FontSize',14,'FontName','Cambria')
        set(cb,'position',[0.9 0.12 0.02 0.76],'linewidth',1.5,'FontSize',12,'FontName','Cambria')
        set(p,'LineWidth',1.5,'Box','on','Ydir','normal','FontSize',12,'FontName','Cambria')
        title(p(2),[name ' TKE Dissipation Rate Timeseries'],'FontSize',14,'FontName','Cambria')
        figname = [name 'TKEtimeseries'];
        export_fig([savefigdir figname],'-jpeg','-nocrop')
    elseif any(isfile) == 0
        continue
    end
end
clear
files = {'F2F2_1TKE.mat';'F2F2_2TKE.mat'};
z = [125 125 80;61 61 61];
datdir1 = 'd:\Projects\Mekong_W2015\DataAnalysis\Paper2\TKE\';
fname = dir(datdir1);fname = {fname.name};

savefigdir = 'd:\Projects\Mekong_W2015\Figures\Paper2\';
for i = 1:length(fname)
    isfile = strcmp(fname{i},files);
    if any(isfile) == 1
        name = regexprep(fname{i},'TKE.mat','');
        load([datdir1 fname{i}])
        
        f1 = figure(i);
        set(f1,'PaperOrientation','portrait',...
            'position',[400 200   1200   500]);
        set(gcf,'color','w','PaperPositionMode','auto')
        c = 1;
        fn = fieldnames(Stat);
        for ii = 1:3
            E = (Stat.(fn{ii}).beam1.E+Stat.(fn{ii}).beam1.E)./2;
            time = Stat.(fn{ii}).time;
            rb =  [0.0400;0.0410;0.0419;0.0429;0.0438;0.0448;0.0457;0.0467;0.0477;0.0486;0.0496;0.0505;...
                0.0515;0.0524;0.0534;0.0543;0.0553;0.0563;0.0572;0.0582;0.0591;0.0601;0.0610;0.0620;...
                0.0630;0.0639;0.0649;0.0658;0.0668;0.0677;0.0687;0.0697;0.0706;0.0716;0.0725];            
            h = (z(c,ii)./1000)-rb;

            %plot
            p(ii) = subplot(1,3,ii);
            imagesc(time,h,E)
            [~,~,~,hr,mi,~] = datevec(time(end)-time(1));
            tint = (hr*60) + mi;                                                    %calculate minutes elapsed between start/stop time of experiment
            tstep = datenum(0,0,0,0,floor(tint/4),0);
            set(gca,'Xtick',time(1):tstep:time(end))
            datetick('x','HH:MM','keepticks','keeplimits')
            if ii == 2
                xlabel(['Time on ' datestr(time(1),'dd/mm/yyyy')],'FontSize',14,'FontName','Cambria')
                set(gca,'Yticklabel',[])
            end
            if ii == 1
                ylabel('Height above bed (m)','FontSize',14,'FontName','Cambria')
            end
            if ii == 3
                set(gca,'Yticklabel',[])
            end
            if c < 2
                c = c+1;
            end
        end
        %global plot adjustments
        set(p(1),'position',[0.08 0.12 0.25 0.76])
        set(p(2),'position',[0.36 0.12 0.25 0.76])
        set(p(3),'position',[0.64 0.12 0.25 0.76])
        cb = colorbar;
        ylabel(cb,'Dissipation Rate (m^2s^-^3)','FontSize',14,'FontName','Cambria')
        set(cb,'position',[0.9 0.12 0.02 0.76],'linewidth',1.5,'FontSize',12,'FontName','Cambria')
        set(p,'LineWidth',1.5,'Box','on','Ydir','normal','FontSize',12,'FontName','Cambria')
        title(p(2),[name ' TKE Dissipation Rate Timeseries'],'FontSize',14,'FontName','Cambria')
        figname = [name 'TKEtimeseries'];
        export_fig([savefigdir figname],'-jpeg','-nocrop')
    elseif any(isfile) == 0
        continue
    end
end