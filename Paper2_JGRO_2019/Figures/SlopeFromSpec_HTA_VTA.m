%Calculate spectra slopes over the time of the deployment for the HTA and
%VTA (SlopeFromSpec version 2). 
clear
close all
datdir = 'd:\Projects\Mekong_W2015\DataAnalysis\Paper2\';
% files = {'HTA_1Vels.mat';'HTA_2Vels.mat';'HTA_4Vels.mat';...
%     'VTA_2vp1Vels.mat';'VTA_2vp2Vels.mat';'VTA_2vp3Vels.mat'};
% Spec = struct();
days = {'day1';'day2';'day3';'day4'};
% for i = 3:6
%     disp(['Loading ' files{i}])
%     load([datdir files{i}])
%     fn = fieldnames(dat);
%     if length(fn) == 3
%         g = 1:3;
%     else
%         g = 1;
%     end
%     if ~isempty(strfind(files{i},'1'))
%         bins = 1:5;
%     else
%         bins = 9:23;
%     end
%     for ii = g
%         time = dat.(fn{ii}).time;    
%         if i < 3
%             heading = 20;
%             ind = 1:length(time)-1;
%             date = days{i};
%         elseif i == 3
%             heading = 20;
%             start = datenum(2015,03,10,15,14,00);
%             stop = datenum(2015,03,10,16,30,00);
%             ind = find(time >= start & time <= stop);
%             time = time(ind);
%             date = days{i};
%         else
%             heading = 96;
%             start = datenum(2015,03,14,07,10,09);
%             stop = datenum(2015,03,14,09,20,09);
%             ind = find(time >= start & time <= stop);
%             time = time(ind);
%             date = days{4};
%         end
%         x = nanmean(dat.(fn{ii}).x(ind,bins),2);
%         y = nanmean(dat.(fn{ii}).y(ind,bins),2);
%         rot = (pi*heading)/180;
%         T = [cos(rot) -sin(rot);...
%             sin(rot) cos(rot)];
%         vels = [x y];
%         V = vels*T';
%         x = V(:,1);y = V(:,2);
%         %loop in time
%         fs = 50;
%         avt = fs*30;
%         nwin = fs*600;
%         swin = fs*10; %10 second averaging window (to smooth)
%         DOF = round((nwin/swin)*2);
%         nsamp = length(time);
%         ind = [1 avt:avt:nsamp];
%         for j = 1:length(ind)
%             if abs(nsamp-ind(j)) < nwin  %skip the last few indexes approaching the end of the t-s
%                 continue
%             else
%                 idx = ind(j):ind(j)+nwin-1;
%             end
%             u = y(idx);
%             v = x(idx);
%             time2 = time(ind(j));
%             
%             [Cuu,F] = pwelch(u,hanning(swin),swin*0.5,nwin,fs); %along shore
%             [Cvv,~] = pwelch(v,hanning(swin),swin*0.5,nwin,fs); %cross shore
%             
%             minf = 4;
%             maxf = 10;
%             fc = find(F >= minf & F <= maxf);
%             %along shore
%             inert = Cuu(fc);ff = F(fc);
%             p = polyfit(log(ff),log(inert),1);
%             Uslope = p(1);
%             %cross shore
%             inert = Cvv(fc);
%             p = polyfit(log(ff),log(inert),1);
%             Vslope = p(1);
% 
% %             figure(1)
% %             y_hat=exp(p(1)*log(ff)+p(2));
% %             loglog(ff,inert,'+','Color','k'), hold on
% %             loglog(ff,y_hat,'--r');
%             Spec.(date).(fn{ii}).Uslope(j) = Uslope;
%             Spec.(date).(fn{ii}).Vslope(j) = Vslope;
%             Spec.(date).(fn{ii}).time(j) = time2;
%         end
%     end
%     clear dat
% end
datdir = 'd:\Projects\Mekong_W2015\DataAnalysis\Paper1\';
% save([datdir 'SpecSlopes'],'Spec','-V7.3')
load([datdir 'SpecSlopes.mat'])
fn = fieldnames(Spec.day1);
for i = 1:4
    for ii = 1:3
        t = Spec.(days{i}).(fn{ii}).time;
        id = find(t<t(1));
        t(id) = NaN;
        bd=isnan(t);
        gd=find(~bd);
        bd([1:(min(gd)-1) (max(gd)+1):end])=0;
        t(bd)=interp1(gd,t(gd),find(bd));
        Spec.(days{i}).(fn{ii}).time = t;
    end
end
%%%Plot Routine
f1 = figure(1);
set(f1,'PaperOrientation','portrait',...
    'position',[400 200   1200   400]);
set(gcf,'color','w','paperpositionmode','auto')
lines = {'--';'-.';':'};
c = [0.7 0.7 0.7;0.5 0.5 0.5;0.1 0.1 0.1];
sp = zeros(1,4);
pp = zeros(1,3);
for i = 1:4
    sp(i) = subplot(1,4,i);
    fn = fieldnames(Spec.(days{i}));
    if i < 4
        symb = {'o';'d';'p'};
        text = {'x = -10 cm';'x = 10 cm';'x = 20 cm'};
    else
        symb = {'o';'s';'d'};
        text = {'z/h_c = 0.04';'z/h_c = 0.6';'z/h_c = 1.25'};
        c = flipud(c);
    end
    for ii = 1:3
        t = Spec.(days{i}).(fn{ii}).time;    
        plot(linspace(t(1),t(end),100),...
        ones(1,100)*-5/3,'-k'),hold on
        n = round(length(t)/10);
        along = smooth(Spec.(days{i}).(fn{ii}).Uslope,n);
        cross = smooth(Spec.(days{i}).(fn{ii}).Vslope,n);
        plot(t,along,lines{ii},...
            'Color',c(ii,:),...
            'LineWidth',1.5),
        markx = t(1:40:end);
        marky = along(1:40:end);
        plot(markx,marky,symb{ii},...
            'markersize',4,...
            'MarkerFaceColor',c(ii,:),...
            'MarkerEdgeColor','k',...
            'LineWidth',1.5)
        plot(t,cross,'-',...
            'Color',c(ii,:),...
            'LineWidth',1.5)
        markx = t(1:40:end);
        marky = cross(1:40:end);
        plot(markx,marky,symb{ii},...
            'markersize',4,...
            'MarkerFaceColor',c(ii,:),...
            'MarkerEdgeColor','k',...
            'LineWidth',1.5)
        pp(ii) = plot(0,0,'-',...
            'marker',symb{ii},...
            'color',c(ii,:),...
            'markersize',4,...
            'MarkerFaceColor',c(ii,:),...
            'MarkerEdgeColor','k',...
            'LineWidth',1.5);
        if i > 2
            tstep = datenum(0,0,0,0,30,0); 
        else
           tstep = datenum(0,0,0,1,0,0); 
        end
        set(gca,'xlim',[t(1) t(end)],...
            'xtick',t(1):tstep:t(end))
        datetick('x','HH:MM','keepticks','keeplimits')
    end
    if i == 3
        leg(1) = legend(pp,text);
    elseif i == 4
        leg(2) = legend(pp,text);        
    end
end
%plot adjustments
set(sp,'ylim',[-2.8 -1])
set(sp(1),'position',[0.09 0.15 0.2 0.7])
set(sp(2),'position',[0.32 0.15 0.2 0.7],...
    'yticklabel',[])
set(sp(3),'position',[0.55 0.15 0.2 0.7],...
    'yticklabel',[])
set(sp(4),'position',[0.78 0.15 0.2 0.7],...
    'yticklabel',[])        
prettyfigures('text',13,'labels',14,'box',1)
set(leg(1),'position',[0.68 0.28 0.02 0.02],'box','off')
set(leg(2),'position',[0.9 0.72 0.02 0.02],'box','off')
xlabel(sp(1),['Time on ' datestr(Spec.day1.vpro1.time(1),...
    'dd-mm')])
xlabel(sp(2),['Time on ' datestr(Spec.day2.vpro1.time(1),...
    'dd-mm')])
xlabel(sp(3),['Time on ' datestr(Spec.day3.vpro1.time(1),...
    'dd-mm')])
xlabel(sp(4),['Time on ' datestr(Spec.day4.vpro1.time(1),...
    'dd-mm')])
ylabel(sp(1),'Spectral Slope (4<f<15 Hz)')
title(sp(1),'HTA, z/h_c = 0.03')
title(sp(2),'HTA, z/h_c = 0.33')
title(sp(3),'HTA, z/h_c = 0.85')
title(sp(4),'VTA')
figdir = 'd:\Projects\Mekong_W2015\Figures\Paper1\';
% export_fig([figdir 'SlopeSpec'],'-png')
