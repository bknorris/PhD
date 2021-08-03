%Calculate spectra slopes over the time of the deployment for the HTA and
%VTA (SlopeFromSpec version 2). 
clear
close all
datdir = 'd:\Mekong_W2015\DataAnalysis\Paper1\';
% files = {'HTA_1Vels.mat';'HTA_2Vels.mat';'HTA_4Vels.mat';...
%     'VTA_2vp1Vels.mat';'VTA_2vp2Vels.mat';'VTA_2vp3Vels.mat'};
% Spec = struct();
% days = {'day1';'day2';'day3';'day4'};
% cutoff = [8.2 8.2 8.2;7 7 7;7 7 7;5.1 0 0 ;5.1 0 0;5.1 0 0]; %from wave cutoff freq (see WaveStatistics_v2)
% for i = 1:6
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
%         disp(fn{ii})
%         time = dat.(fn{ii}).time;    
%         if i == 1
%             heading = 20;
%             ind = 1:length(time)-1;
%             date = days{i};
%         elseif i == 2
%             heading = 20;
%             start = datenum(2015,03,08,13,48,00);
%             stop = datenum(2015,03,08,16,48,00);
%             ind = find(time >= start & time <= stop);
%             time = time(ind);
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
%         z = nanmean((dat.(fn{ii}).z1(ind,bins)+dat.(fn{ii}).z1(ind,bins))./2,2);
%         rot = (pi*heading)/180;
%         T = [cos(rot) -sin(rot);...
%             sin(rot) cos(rot)];
%         vels = [x y];
%         V = vels*T';
%         x = V(:,1);y = V(:,2);
%         %loop in time
%         fs = 50;
%         avt = fs*60;
%         nwin = fs*300;
%         swin = fs*20; %10 second averaging window (to smooth)
%         DOF = round((nwin/swin)*2);
%         nsamp = length(time);
%         ind = [1 avt:avt:nsamp];
% %         c = jet(length(ind));
%         for j = 1:length(ind)
%             if abs(nsamp-ind(j)) < nwin  %skip the last few indexes approaching the end of the t-s
%                 continue
%             else
%                 idx = ind(j):ind(j)+nwin-1;
%             end
%             u = y(idx);
%             v = x(idx);
%             w = z(idx);
%             time2 = time(ind(j));
%             
%             [Cuu,F] = pwelch(u,hanning(swin),swin*0.5,nwin,fs); %along shore
%             [Cvv,~] = pwelch(v,hanning(swin),swin*0.5,nwin,fs); %cross shore
%             [Cww,~] = pwelch(w,hanning(swin),swin*0.5,nwin,fs); %cross shore
% %             figure(1)
% %             loglog(F,Cvv,'color',c(j,:)),hold on
%             
%             minf = cutoff(i,ii);
%             maxf = 15;
%             if j == 1
%                 disp(['Applying cutoff: ' sprintf('%0.1f',minf) ' < F < ' sprintf('%0.1f',maxf)])
%             end
%             fc = find(F >= minf & F <= maxf);
%             %along shore
%             inert = Cuu(fc);ff = F(fc);
%             p = polyfit(log(ff),log(inert),1);
%             Uslope = p(1);
%             %cross shore
%             inert = Cvv(fc);
%             p = polyfit(log(ff),log(inert),1);
%             Vslope = p(1);
%             %vertical
%             inert = Cww(fc);
%             p = polyfit(log(ff),log(inert),1);
%             Wslope = p(1);
% 
% %             figure(2)
% %             y_hat=exp(p(1)*log(ff)+p(2));
% %             loglog(ff,inert,'+','Color','k'), hold on
% %             loglog(ff,y_hat,'--r');
%             
%             Spec.(date).(fn{ii}).Uslope(j) = Uslope;
%             Spec.(date).(fn{ii}).Vslope(j) = Vslope;
%             Spec.(date).(fn{ii}).Wslope(j) = Wslope;
%             Spec.(date).(fn{ii}).time(j) = time2;
%         end
%     end
%     clear dat
% end

datdir = 'd:\Mekong_W2015\DataAnalysis\Paper2\';
% save([datdir 'SpecSlopes'],'Spec','-V7.3')

load([datdir 'SpecSlopes (2).mat'])
% fn = fieldnames(Spec.day1);
% for i = 1:4
%     for ii = 1:3
%         t = Spec.(days{i}).(fn{ii}).time;
%         id = find(t<t(1));
%         t(id) = NaN;
%         bd=isnan(t);
%         gd=find(~bd);
%         bd([1:(min(gd)-1) (max(gd)+1):end])=0;
%         t(bd)=interp1(gd,t(gd),find(bd));
%         Spec.(days{i}).(fn{ii}).time = t;
%     end
% end

%%%Plot Routine
f1 = figure(1);
set(f1,'PaperOrientation','portrait',...
    'position',[400 100   1000   400]);
set(gcf,'color','w','paperpositionmode','auto')   
symb = {'p';'d';'^'};
text = {'x = 20 cm';'x = 10 cm';'z/h_c = 0.04'};
days = {'day1';'day2';'day4'};
vps = [3 2 1];
c = [0.2 0.2 0.2;0.5 0.5 0.5;0 0 0];
sp = zeros(1,3);
leg = zeros(1,3);
for i = 1:3
    sp(i) = subplot(1,3,i);
    fn = fieldnames(Spec.(days{i}));
    t = Spec.(days{i}).(fn{vps(i)}).time;
    t(t < 1E5) = NaN;t = fixgaps(t);
    plot(linspace(t(1),t(end),100),...
        ones(1,100)*-5/3,'-k'),hold on
    n = round(length(t)/10);
    along = fastsmooth(Spec.(days{i}).(fn{vps(i)}).Uslope,n,1,1);
    cross = fastsmooth(Spec.(days{i}).(fn{vps(i)}).Vslope,n,1,1);
    vert = fastsmooth(Spec.(days{i}).(fn{vps(i)}).Wslope,n,1,1);
    plot(t,along,'--',...
        'Color',c(i,:),...
        'LineWidth',1.5),
    markx = t(1:30:end);
    marky = along(1:30:end);
    plot(markx,marky,symb{i},...
        'markersize',6,...
        'MarkerFaceColor',c(i,:),...
        'MarkerEdgeColor','k',...
        'LineWidth',1.5)
    plot(t,cross,'-',...
        'Color',c(i,:),...
        'LineWidth',1.5)
    markx = t(1:30:end);
    marky = cross(1:30:end);
    plot(markx,marky,symb{i},...
        'markersize',6,...
        'MarkerFaceColor',c(i,:),...
        'MarkerEdgeColor','k',...
        'LineWidth',1.5)
    plot(t,vert,':',...
        'Color',c(i,:),...
        'LineWidth',1.5)
    markx = t(1:30:end);
    marky = vert(1:30:end);
    plot(markx,marky,symb{i},...
        'markersize',6,...
        'MarkerFaceColor',c(i,:),...
        'MarkerEdgeColor','k',...
        'LineWidth',1.5)
    pp = plot(0,0,'-',...
        'marker',symb{i},...
        'color',c(i,:),...
        'markersize',6,...
        'MarkerFaceColor',c(i,:),...
        'MarkerEdgeColor','k',...
        'LineWidth',1.5);
        if i > 1
            tstep = datenum(0,0,0,0,30,0); 
        else
           tstep = datenum(0,0,0,1,0,0); 
        end
        if i < 3
            set(gca,'xlim',[t(20) t(end-20)],...
                'xtick',t(20):tstep:t(end-20))
        else
            set(gca,'xlim',[t(1) t(end)],...
                'xtick',t(1):tstep:t(end))
        end
        datetick('x','HH:MM','keepticks','keeplimits')
        leg(i) = legend(pp,text{i});
end
%plot adjustments
set(sp,'ylim',[-3.4 -1])
set(sp(1),'position',[0.1 0.15 0.26 0.7])
set(sp(2),'position',[0.4 0.15 0.26 0.7],...
    'yticklabel',[])
set(sp(3),'position',[0.72 0.15 0.26 0.7],...
    'yticklabel',[])

prettyfigures('text',13,'labels',14,'box',1)
set(leg,'box','off')
xlabel(sp(1),['Time on ' datestr(Spec.day1.vpro1.time(1),...
    'dd-mm')])
xlabel(sp(2),['Time on ' datestr(Spec.day2.vpro1.time(1),...
    'dd-mm')])
xlabel(sp(3),['Time on ' datestr(Spec.day4.vpro1.time(1),...
    'dd-mm')])
ylabel(sp(1),'Spectral Slope (f_c < f < 15 Hz)')
title(sp(1),'HTA1, z/h_c = 0.03')
title(sp(2),'HTA2, z/h_c = 0.33')
title(sp(3),'VTA')
figdir = 'c:\Users\Bnorr\Documents\GradSchool\Writing\JGR_2018_Norris\Figures\Draft2_Figures\Versions\V2\';
export_fig([figdir 'SlopeSpec_V3'],'-pdf')
