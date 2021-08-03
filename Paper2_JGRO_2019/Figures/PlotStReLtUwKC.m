%Plot St, Re, Uw and Uc, Lt and KC for the HTA and VTA experiments
clear
load('D:\Projects\Mekong_W2015\DataAnalysis\Paper2\WVStats_HTA_VTA.mat')
wv = data;
load('D:\Projects\Mekong_W2015\DataAnalysis\Paper2\HTAVTA_UwUcStLt.mat')
fn = fieldnames(data);sn = fieldnames(wv);
c = 1;
for i = 1:length(fn)
    dn = fieldnames(data.(fn{i}));
    disp(fn{i})
    for ii = 1:3
        f1 = figure(c);
        set(f1,'PaperOrientation','portrait',...
            'position',[500 100   800   1000]);
        suptitle([upper(fn{i}) ' ' upper(dn{ii})])
        time = data.(fn{i}).(dn{ii}).time;time(time < 1) = NaN;
        subplot(511)
        pp(1) = plot(time,data.(fn{i}).(dn{ii}).Uw,'-r');hold on
        pp(2) = plot(time,data.(fn{i}).(dn{ii}).Uc,'-b');
        legend(pp,{'U_w','U_c'})
        set(gca,'xticklabel',[])
        ylabel('m/s')
        title('Wave Velocities')
        subplot(512)
        plot(time,data.(fn{i}).(dn{ii}).Re,'-k');
        title('Reynolds Number')
        set(gca,'xticklabel',[])
        ylabel('Re')
        subplot(513)
        plot(time,data.(fn{i}).(dn{ii}).St(:,1),'-r');hold on
        plot(time,data.(fn{i}).(dn{ii}).St(:,2),'-b');
        plot(time,data.(fn{i}).(dn{ii}).St(:,3),'-g');
        plot(time,0.2*ones(length(time),1),'-k');
        title('Strouhal Number')
        set(gca,'xticklabel',[])
        ylabel('St')
        subplot(514)
        plot(time,data.(fn{i}).(dn{ii}).Lt(:,1),'-r');hold on
        plot(time,data.(fn{i}).(dn{ii}).Lt(:,2),'-b');
        plot(time,data.(fn{i}).(dn{ii}).Lt(:,3),'-g');
        plot(time,data.(fn{i}).(dn{ii}).d,'-k');
        title('Integral Length Scale')
        set(gca,'xticklabel',[])
        ylabel('Lt (m)')
        %calculate KC = Uw*T/d
        %spline wave period to Uw
        Tr = wv.(sn{i}).Tr;Tr(isnan(Tr)) = 0;
        T = interp1(1:length(wv.(sn{i}).time2),Tr,1:length(time));
        d = data.(fn{i}).(dn{ii}).d;
        KC = (data.(fn{i}).(dn{ii}).Uw.*T)./d;
        subplot(515)
        plot(time,KC,'-r');hold on
        plot(time,6*ones(length(time),1),'-k');
        title('Keulegan-Carpenter Number')
        datetickzoom('x','HH:MM','keepticks','keeplimits')
        xlabel(['Time on ' datestr(time(1),'dd-mm-yyyy')])
        ylabel('KC')
        c = c+1;
        
        %Disp min, mean and max of parameters
        if i == 1
            sta = 128;
            stp = length(KC);
        elseif i == 2;
            sta = 1;
            stp = 1137;
        elseif i == 4;
            sta = 145;
            stp = 971;
        else
            sta = 1;
            stp = length(KC);
        end
        Uwmi = nanmin(data.(fn{i}).(dn{ii}).Uw(sta:stp));
        Uwm = nanmean(data.(fn{i}).(dn{ii}).Uw(sta:stp));
        Uwma = nanmax(data.(fn{i}).(dn{ii}).Uw(sta:stp));
        
        Ucmi = nanmin(data.(fn{i}).(dn{ii}).Uc(sta:stp));
        Ucm = nanmean(data.(fn{i}).(dn{ii}).Uc(sta:stp));
        Ucma = nanmax(data.(fn{i}).(dn{ii}).Uc(sta:stp));
        
        Ltmi = nanmin(data.(fn{i}).(dn{ii}).Lt(sta:stp));
        Ltm = nanmean(data.(fn{i}).(dn{ii}).Lt(sta:stp));
        Ltma = nanmax(data.(fn{i}).(dn{ii}).Lt(sta:stp));
        
        Remi = nanmin(data.(fn{i}).(dn{ii}).Re(sta:stp));
        Rem = nanmean(data.(fn{i}).(dn{ii}).Re(sta:stp));
        Rema = nanmax(data.(fn{i}).(dn{ii}).Re(sta:stp));
        
        Stmi = nanmin(data.(fn{i}).(dn{ii}).St(sta:stp));
        Stm = nanmean(data.(fn{i}).(dn{ii}).St(sta:stp));
        Stma = nanmax(data.(fn{i}).(dn{ii}).St(sta:stp));
        
        disp(dn{ii})
        disp(['Uw min: ' num2str(Uwmi)])
        disp(['Uw mean: ' num2str(Uwm)])
        disp(['Uw max: ' num2str(Uwma)])
        
        disp(['Uc min: ' num2str(Ucmi)])
        disp(['Uc mean: ' num2str(Ucm)])
        disp(['Uc max: ' num2str(Ucma)])
        
        disp(['Lt min: ' num2str(Ltmi)])
        disp(['Lt mean: ' num2str(Ltm)])
        disp(['Lt max: ' num2str(Ltma)])
        
        disp(['Re min: ' num2str(Remi)])
        disp(['Re mean: ' num2str(Rem)])
        disp(['Re max: ' num2str(Rema)])
        
        disp(['St min: ' num2str(Stmi)])
        disp(['St mean: ' num2str(Stm)])
        disp(['St max: ' num2str(Stma)])
        fprintf('\n')
    end
end


