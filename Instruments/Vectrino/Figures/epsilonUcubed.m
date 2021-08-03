%Does epsilon scale with u^3? - a test of... uh...

clear
load('C:\Users\bkn5\Projects\Mekong_W2015\DataAnalysis\Spectra\Paper1\QCd\HTAday1Vels.mat')
statdir = 'C:\Users\bkn5\Projects\Mekong_W2015\DataAnalysis\WavePowerFlux\FSS3Turbulence\';
figdir = 'C:\Users\bkn5\Projects\Mekong_W2015\Figures\Turbulence\';
fn = fieldnames(HTA);    
bin = 13; %bin # to plot
tstr = 'Near Bottom Vectrinos '; 
%title string text

count = 1;
f = figure;
set(f,'PaperOrientation','portrait',...
    'position',[400 200   600   500]);
set(gcf,'color','w','PaperPositionMode','auto')
markers = {'^','d','o','s','x'};

c = linspecer(5);
for i = 2:4
    files = dir([statdir '*070315*']);
    fname = {files.name};
    load([statdir fname{count}])
    
    x = HTA.(fn{i}).x;
    y = HTA.(fn{i}).y;
    intv = 10;
    fs = 50;
    avt =  60*fs*intv; %samples/interval for averaging
    ind = [1 avt:avt:length(x)];
    for j = 1:length(ind)-1
        u = x(ind(j):ind(j+1),bin); %bin 3
        u = abs(u);
        ucb = u.^3;
        U(:,j) = mean(ucb);
    end
    %take epsilon to be the average between the TKE estimates of beams 1 and
    %beam 3: the cross-shore beams
    epsilon = (STAT.TKE.Beam1(bin,:)+STAT.TKE.Beam3(bin,:))./2;
    if count == 3 %last instrument
        %deposition beneath the instrument is making the TKE estimates
        %higher than others; eliminate outliers
        epsilon = epsilon(1:12);
        U = U(1:12);
    end
    pf = polyfit(U,epsilon,1);
    pv = polyval(pf,U);
    %calculate residuals as a measure of fit
    resid = epsilon-pv;
    SSresid = sum(resid.^2);
    SStotal = (length(epsilon)-1)*var(epsilon);
    rsq = 1-SSresid/SStotal;
    disp(['R-squared value: ' num2str(rsq)])
    
    g(count) = plot(U,epsilon,markers{count},'Color',c(count,:),'MarkerSize',8); hold on
    xb = linspace(min(U),6E-3,length(pv));
    yb = pf(1)*xb+pf(2);
    plot(xb,yb,'-','LineWidth',1,'Color',c(count,:))
    count = count+1;
end

load('C:\Users\bkn5\Projects\Mekong_W2015\DataAnalysis\Spectra\Paper1\QCd\F2F2day2vels.mat')
load('c:\Users\bkn5\Projects\Mekong_W2015\DataAnalysis\WavePowerFlux\FSS3Turbulence\Stat_VP1_060315.mat')
x = F2F2.vpro1.x;
y = F2F2.vpro1.y;
intv = 10;
fs = 50;
avt =  60*fs*intv; %samples/interval for averaging
ind = [1 avt:avt:length(x)];
for j = 1:length(ind)-1
    u = x(ind(j):ind(j+1),bin); %bin 3
    u = abs(u);
    ucb = u.^3;
    U(:,j) = mean(ucb);
end
%take epsilon to be the average between the TKE estimates of beams 1 and
%beam 3: the cross-shore beams
epsilon = (STAT.TKE.Beam1(bin,:)+STAT.TKE.Beam3(bin,:))./2;U = U(1:20);
pf = polyfit(U,epsilon,1);
pv = polyval(pf,U);
%calculate residuals as a measure of fit
resid = epsilon-pv;
SSresid = sum(resid.^2);
SStotal = (length(epsilon)-1)*var(epsilon);
rsq = 1-SSresid/SStotal;
disp(['R-squared value: ' num2str(rsq)])

g(count) = plot(U,epsilon,markers{count},'Color',c(count,:),'MarkerSize',8); hold on
xb = linspace(min(U),6E-3,length(pv));
yb = pf(1)*xb+pf(2);
plot(xb,yb,'-','LineWidth',1,'Color',c(count,:))

count = count+1;

load('C:\Users\bkn5\Projects\Mekong_W2015\DataAnalysis\Spectra\Paper1\VTAvelocities.mat')
load('c:\Users\bkn5\Projects\Mekong_W2015\DataAnalysis\WavePowerFlux\DPS2Turbulence\Stat_VP1_140315.mat')
x = VTA.vpro1.x;
y = VTA.vpro1.y;
intv = 10;
fs = 50;
avt =  60*fs*intv; %samples/interval for averaging
ind = [1 avt:avt:length(x)];
for j = 1:length(ind)-1
    u = x(ind(j):ind(j+1),bin); %bin 3
    u = abs(u);
    ucb = u.^3;
    U(:,j) = mean(ucb);
end
%take epsilon to be the average between the TKE estimates of beams 1 and
%beam 3: the cross-shore beams
epsilon = (STAT.TKE.Beam1(bin,:)+STAT.TKE.Beam3(bin,:))./2;epsilon = epsilon(3:end);U = U(1:12);
pf = polyfit(U,epsilon,1);
pv = polyval(pf,U);
%calculate residuals as a measure of fit
resid = epsilon-pv;
SSresid = sum(resid.^2);
SStotal = (length(epsilon)-1)*var(epsilon);
rsq = 1-SSresid/SStotal;
disp(['R-squared value: ' num2str(rsq)])

g(count) = plot(U,epsilon,markers{count},'Color',c(count,:),'MarkerSize',8); hold on
xb = linspace(min(U),6E-3,length(pv));
yb = pf(1)*xb+pf(2);
plot(xb,yb,'-','LineWidth',1,'Color',c(count,:))

leg = legend(g,{'HTA VP1';'HTA VP2';'HTA VP3';'F2F2 VP1';'VTA VP1'});
xlabel('|U^3|')
ylabel('\epsilon (W/kg)')
title([tstr])
figname = ([regexprep(tstr,' ','') 'UcubedE']);
export_fig([figdir figname],'-jpeg')

