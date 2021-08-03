%Plot turbulence for the HTA & VTA experiments as a longitudinal profile across
%the canopy (12/12/16). Epsilon values are normalized by the wave energy
%flux
%NOTE: change "epsilon_t/epsilon_w" to "\widetilde{\epsilon}" if further changes need to be made.
clear
close all
figdir = 'c:\Users\Bnorr\Documents\GradSchool\Papers\JGR_2018_Norris\Figures\Draft2_Figures\Versions\V2\';
ddir = 'd:\Mekong_W2015\DataAnalysis\Paper2\TKE\Vertical\';
wdir = 'd:\Mekong_W2015\DataAnalysis\Paper2\WPF\';
dn = {'day1';'day2';'day3';'day4'};
files = {'HTA_1TKE_bdadj';'HTA_2TKE';'HTA_4TKE';'VTA_2TKE_bdadj'};
wfiles = {'WpfV5108_HTA_1';'Wpf5116_HTA_1';...
    'WpfHR3_HTA_2';'Wpf5116_HTA_2';...
    'WpfHR3_HTA_4';'Wpf5116_HTA_4';...
    'WpfV5109_VTA';'Wpf5116_VTA'};
wf = [1 2;3 4;5 6;7 8];
delx = [65 56 56 31]; %m
theta = [-54 0 0 0]; %deg
rho = 1011.4; %kg/m^3

for i = 1:4
    load([ddir files{i} '.mat'])
    w1 = load([wdir wfiles{wf(i,1)}]);
    w2 = load([wdir wfiles{wf(i,2)}]);
    fn = fieldnames(Stat);
    for ii = 1:3
        %Average data into 10-minute chunks
        n = length(Stat.(fn{ii}).time);
        fs = 0.1; %Hz, epsilon calculated every 10 sec
        avt = 10*60*fs;
        idx = [1 avt:avt:n];
        etmean = zeros(30,length(idx)-1);
        edmean = zeros(length(idx)-1,1);
        enorm = zeros(length(idx)-1,1);
        time = zeros(length(idx)-1,1);
        h = w2.WPF.h;
        wpf = abs((w1.WPF.F-w2.WPF.F)./(delx(i).*cosd(theta(i))));
        for j = 1:length(idx)-1
            %use top 5 bins for HTA day 1, otherwise use middle 15 bins.
            %VTA uses whole profile
            if i == 1
                bins = 1:5;
            elseif i < 4
                bins = 9:23;
            else
                bins = 1:30;
            end
            ids = idx(j):idx(j+1);
            dmean = (nanmean(Stat.(fn{ii}).z1.E(bins,ids))+...
                nanmean(Stat.(fn{ii}).z2.E(bins,ids)))./2; %depth mean
            tmean = (nanmean(Stat.(fn{ii}).z1.E(bins,ids),2)+...
                nanmean(Stat.(fn{ii}).z2.E(bins,ids),2))./2; %time mean
            edmean(j) = mean(dmean);
            etmean(bins,j) = (tmean.*mean(h(ids))*rho)./mean(wpf(ids));
            enorm(j) = mean((dmean'.*h(ids)*rho)./wpf(ids));
            time(j) = Stat.(fn{ii}).time(idx(j));
        end
        data.(dn{i}).(fn{ii}).time = time;
        data.(dn{i}).(fn{ii}).edmean = edmean;
        data.(dn{i}).(fn{ii}).etmean = etmean;
        data.(dn{i}).(fn{ii}).enorm = enorm;
    end
    if i < 4
        v1 = data.(dn{i}).vpro1.enorm;
        v2 = data.(dn{i}).vpro2.enorm;
        v3 = data.(dn{i}).vpro3.enorm;
        array = [v1 v2 v3];
        array(:,2) = (v1+v2)./2;
        data.(dn{i}).array = array;
    end
end
dist = [-10 10 20];

%%%Plot Routine
f1 = figure(1);
set(f1,'PaperOrientation','portrait',...
    'position',[400 100   1000   700]);
symb = {'o';'d';'p'};
cl = flipud([0.1 0.1 0.1;0.5 0.5 0.5;0.7 0.7 0.7]);
pl = [3 2 1];
txt = {'z/h_c = 0.03';'z/h_c = 0.33';'z/h_c = 0.85'};
%set colors based on timing of the tide, between low tide and the end of
%the longest VP deployment
t1 = datenum(2015,03,08,13,20,00);
tend = datenum(2015,03,08,19,00,00);
dt = abs(t1-tend);[~,~,~,h,mn,~] = datevec(dt);
telap = h*60+mn;
time = [1 10:10:telap];
toff = [0 20 40];
c = brewermap(length(time),'PuBu');
sp = zeros(3,1);
%Plot HTA
for i = 1:3
    sp(i) = subplot(6,1,pl(i));
    n = length(data.(dn{i}).array);
    [~,~,~,h,mn,s] = datevec(data.(dn{i}).vpro1.time(1));
    tstart = h*60+mn+(s/60);
    [~,~,~,h,mn,s] = datevec(t1);
    tidest = h*60+mn+(s/60)+toff(i);
    dt = tstart-tidest;
    tmp = abs(time-dt);
    [~,id] = min(tmp);
    cc = c(id:id+n,:);
    for ii = 1:n
        plot(dist,log10(data.(dn{i}).array(ii,:)),...
            'color',cc(ii,:),...
            'linewidth',1.5),hold on
        for j = 1:3
            plot(dist(j),log10(data.(dn{i}).array(ii,j)),symb{j},...
                'markeredgecolor','k',...
                'markerfacecolor',cl(j,:),...
                'linewidth',1.5)
        end
    end
    plot(zeros(5,1),linspace(-4,1,5),'--k',...
        'linewidth',1.5), hold off
    if i < 3
        text(-10,0.24,txt{i})
    else
        text(-10,0.04,txt{i})
    end
    grid on
    set(gca,'gridlinestyle',':')
end
%Plot VTA
pl = [6 5 4];
sp2 = zeros(3,1);
rb = (1:30)./1000;
vph = [0.07 0.416 0.806];
hc = 0.6;
for i = 1:3
    sp2(i) = subplot(6,1,pl(i));
    n = length(data.day4.(fn{i}).time);
    [~,~,~,h,mn,s] = datevec(data.day4.(fn{i}).time(1));
    tstart = h*60+mn+(s/60);
    [~,~,~,h,mn,s] = datevec(t1);
    tidest = h*60+mn+(s/60)-540;
    dt = tstart-tidest;
    tmp = abs(time-dt);
    [~,id] = min(tmp);
    cc = c(id:id+n,:);
    for ii = 1:n
        etmean = log10(data.day4.(fn{i}).etmean(:,ii));
        yy = (vph(i)-0.04-rb)./hc;
        plot(etmean,yy,'-',...
            'Color',cc(ii,:),...
            'LineWidth',1.5), hold on
        markx = etmean(1:4:end);
        marky = yy(1:4:end);
        plot(markx,marky,'o',...
            'linewidth',1.5,...
            'markersize',6,...
            'markerfacecolor','w',...
            'markeredgecolor',cc(ii,:))
    end
    grid on
    set(gca,'gridlinestyle',':')
end    
caxis([time(1) time(end)])
colormap(c);
cb = colorbar;
set(sp,'xtick',-10:5:20,...
    'xlim',[-11 21],...
    'ylim',[-4 1],...
    'ytick',-4:1:1)
% set([sp(1) sp(2)],'xtick',-10:5:20,...
%     'xlim',[-11 21],...
%     'ylim',[-0.02 0.3],...
%     'ytick',0:0.1:0.3)
% set(sp(3),'xtick',-10:5:20,...
%     'xlim',[-11 21],...
%     'ylim',[-0.0025 0.05],...
%     'ytick',0:0.025:0.05)
set(sp2,'xlim',[-4 0])
set(sp2(1),'ylim',[0 0.05],'ytick',0:0.02:0.05)
set(sp2(2),'ylim',[0.58 0.63],'ytick',0.59:0.02:0.65)
set(sp2(3),'ylim',[1.23 1.28],'ytick',1.22:0.02:1.28)
set([sp(2) sp(3)],'xticklabel',[])
set(sp(1),'position',[0.1 0.1 0.38 0.26])
set(sp(2),'position',[0.1 0.39 0.38 0.26])
set(sp(3),'position',[0.1 0.68 0.38 0.26])
set([sp2(2) sp2(3)],'xticklabel',[])
set(sp2(1),'position',[0.58 0.1 0.25 0.26])
set(sp2(2),'position',[0.58 0.39 0.25 0.26])
set(sp2(3),'position',[0.58 0.68 0.25 0.26])
set(cb,'position',[0.86 0.1 0.02 0.84],...
    'ytick',0:50:time(end))
ylabel(cb,'Time elapsed after low tide (min)')
ylabel(sp(2),'log_1_0(\epsilon_T/\epsilon_w)')
xlabel(sp(1),'Cross shore distance (cm)')
ylabel(sp2(2),'z/h_c')
xlabel(sp2(1),'log_1_0(\epsilon_T/\epsilon_w)')
title(sp(3),'HTA'),title(sp2(3),'VTA')
prettyfigures('text',13,'labels',14,'box',1)
export_fig([figdir 'HTA_VTATurbProfile_v3'],'-pdf','-nocrop')