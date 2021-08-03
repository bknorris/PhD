%Calculate the variance of the bottom trace; plot U (interval averaged
%velocity), bottom trace, and var(bottom trace) for all three VPs, per
%experiment. See Staudt et al. 2017 for reference. Also plot the mean of
%the Amplitude (backscatter) versus bottom trace
%
% Paper 3: Sediment Motion in Mangroves
%
% This is Version 1.0 of this script
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear
close all

expname = 'F2F_2';
files = {'F2F_2_bdtrace.mat';'F2F_2_Vels.mat';'F2F_2_Sen.mat'}; 
hoffset = [20 20 20]; %deg; offset between mag N and transect for rotations
bins = 13:17; %bins for averaging velocity
start = datenum('29-Sep-2014 00:00:00');
stop = datenum('31-Sep-2014 00:00:00');
%%%
dir1 = 'd:\Projects\Mekong_W2015\DataAnalysis\Paper3\BottomTrackCPSD\';
fn1 = whos('-file',[dir1 files{1}]);fn1 = {fn1.name};
trace = matfile([dir1 files{1}]);
%%%
dir2 = 'd:\Projects\Mekong_F2014\DataAnalysis\Paper3\VPs\';
vels = matfile([dir2 files{2}]);sen = matfile([dir2 files{3}]);
if strfind(expname,'VTA')
    ni = 1;
else
    ni = 1:3;
end
%Rangebins height
rb = [0.0400,0.0410,0.0419,0.0429,0.0439,...
0.0448,0.0458,0.0467,0.0477,0.0487,...
0.0496,0.0506,0.0516,0.0525,0.0535,...
0.0545,0.0554,0.0564,0.0574,0.0583,...
0.0593,0.0602,0.0612,0.0622,0.0631,...
0.0641,0.0651,0.0660,0.0670,0.0680,...
0.0689,0.0699,0.0708,0.0718,0.0728];
sdir = 'd:\Projects\Mekong_W2015\Figures\Paper3\BottomTracking\';

%Figure 1: Amplitude w/ Bottom trace
for i = ni
    bd = trace.(fn1{i}); %bottom trace
    cp1 = find(bd.time >= start & bd.time <= stop);
    s = sen.(fn1{i}); %sen file
    time = s.time;
    cp2 = find(time >= start & time <= stop);
    time = time(cp2);
    nsamp = length(time);
    amp = (s.Amp1+s.Amp2+s.Amp3+s.Amp4)./4;
    amp = amp(cp2,:);

    %Interp bdh to 50Hz
    bdh = interp1(bd.time(cp1),bd.bdist(cp1),time);
    bid = find(~isnan(bdh),1,'first');bdi = bdh(bid);
    bdh = bdi-bdh;
    clear s %save memory
    
    %%%
    
    f1 = figure(i);
    set(f1,'PaperOrientation','portrait',...
        'position',[400 100   800 500]);
	z = bdi-rb;
    imagesc(time,z,amp'),hold on
    caxis([-25 -5])
    cb = colorbar;
    set(gca,'ydir','normal')
    
    
    %Find max(amp) near bdh
    amax = zeros(nsamp,1);
    for ii = 1:nsamp
        id = find(z<=bdh(ii)+0.002 & z>=bdh(ii)-0.002);
        if isempty(id)
            continue
        end
        [~,aid] = max(amp(ii,id));
        amax(ii) = z(id(aid));
    end
    am = decimate(amax,5); %downsample to 10Hz
    am = my_running_median(am,512);
    lt = length(am)-length(bd.time(cp1));
    am = interp1(bd.time(cp1),am(1:end-lt),time); %interp to 50Hz
    plot(time,am,'k',...
        'linewidth',1.5)    
    plot(time,bdh,'w',...
        'LineWidth',1.5)
    title([expname ' - ' upper(fn1{i}) ': Amplitude & Bed Trace'])
    ylabel('Height above bed (m)')
    xlabel(['Time on ' datestr(time(1),'dd-mm-yy')])
    datetick('x','HH:MM','keepticks','keeplimits')
    ylabel(cb,'Amplitude (dB)')
    clear am amax amp id aid
end 

%Figure (2): Bottom Trace, U and Var(Bottom Trace)
sp = zeros(3,1);pl = zeros(3,1);
f2 = figure(4);
set(f2,'PaperOrientation','portrait',...
    'position',[400 100   500  1000]);
hold on
if max(size(ni)) == 3
    cl = brewermap(5,'Blues');
    cl = cl(3:5,:);
    csp = [0.6 0.6 0.6;...
        0.4 0.4 0.4;...
        0.1 0.1 0.1];
else
    cl = [0 0 0];
    csp = [1 0 0];
end
miU = zeros(max(size(ni)),1);
maU = zeros(max(size(ni)),1); %scaling for U
for i = ni
    %Calculate Var(Bottom Trace)
    bd = trace.(fn1{i}); %bottom trace
    cp1 = find(bd.time >= start & bd.time <= stop);
    v = vels.(fn1{i});
    time = v.time;
    cp2 = find(time >= start & time <= stop);
    time = time(cp2);
    x = nanmean(v.x(cp2,bins),2);x(x==0)=NaN;
    y = nanmean(v.y(cp2,bins),2);y(y==0)=NaN;
    
    %Interp bdh to 50Hz
    bdh = interp1(bd.time(cp1),bd.bdist(cp1),time);
    bid = find(~isnan(bdh),1,'first');
    bdh = bdh(bid)-bdh;
    clear v; %save memory  

    sp(1) = subplot(3,1,1);
    %Bottom Trace
    plot(time,bdh,'color',cl(i,:),'linewidth',1.5),hold on
    t1 = time(1);t2 = time(end);

    window = 180; %3 min averaging interval
    step = 10; %10 second step
    fs = 50;
    avt = step*fs; %samples/step
    nwin = window*fs; %samples/window
    ind = [1 avt:avt:nsamp];
    %Initialize Variables
    s_bint = zeros(length(ind),1);
    U = zeros(length(ind),1);
    sT = zeros(length(ind),1);
    for ii = 1:length(ind)                                                  
        if abs(nsamp-ind(ii)) < nwin
            sT(ii) = time(idx(1));
            continue
        else
            idx = ind(ii):ind(ii)+nwin-1;
            X = x(idx);
            Y = y(idx);
            dbi = bdh(idx);dbi(isnan(dbi)) = [];
        end
        %Rotate to cross-shore and along-shore
        rot = hoffset(i)*pi/180;
        X = X.*(ones(size(X))*cos(rot)) + ...
            Y.*(ones(size(Y))*sin(rot));
        Y = -Y.*(ones(size(Y))*sin(rot)) + ...
             X.*(ones(size(X))*cos(rot));
        %Calculate sigma^2_b,int (e.g. Staudt et al. 2017)
        s_bint(ii) = var(dbi);
        U(ii) = mean(sqrt(X.^2+Y.^2));
        sT(ii) = time(idx(1));
    end
    s_bint(s_bint==0) = NaN;
    U(U==0) = NaN;
    s_bnorm = s_bint./(U.*window);
    %Average s_bnorm
    avt = 90; %15 min interval where step = 10s
    lt = length(s_bnorm);
    s_bnave = zeros(floor(lt/avt),1);
    sTb = zeros(floor(lt/avt),1);
    for ii = 1:floor(lt/avt)
        s_bnave(ii) = mean(s_bnorm((ii-1)*avt+1:ii*avt));
        sTb(ii) = sT((ii-1)*avt+1);
    end
    sp(2) = subplot(3,1,2);
    %Bottom Trace
    if i == 1
        plot(sT,zeros(length(sT),1),'--',...
            'LineWidth',1.5,'Color',[0.5 0.5 0.5]),hold on
    end
    plot(sT,U,'color',cl(i,:),'linewidth',1.5)
    miU(i) = min(U)-0.01;maU(i) = max(U)+0.01;
    
    sp(3) = subplot(3,1,3);
    pl(i) = plot(sT,s_bnorm,'color',cl(i,:),'linewidth',1.5);hold on      
    plot(sTb,s_bnave,'--','marker','o','markersize',6,...
        'color',csp(i,:),'linewidth',1.5) 
end

%Global adjustments
set([sp(1) sp(2)],'xticklabel',[])
set(sp,'xlim',[t1 t2])
set(sp(2),'ylim',[min(miU) max(maU)])
set(sp(3),'yscale','log')
datetick('x','HH:MM','keepticks','keeplimits')
xlabel(sp(3),['Time on ' datestr(time(1),'dd-mm-yy')])
ylabel(sp(3),'$\frac{\sigma_{b,int}^{2}}{U\Delta t}$',...
    'interpreter','latex')
ylabel(sp(2),'$\sqrt{u^{2}+v^{2}}$   (m/s)',...
    'interpreter','latex')
ylabel(sp(1),'Bed Elevation (m)')
set(sp(3),'position',[0.16 0.1 0.8 0.25])
set(sp(2),'position',[0.16 0.4 0.8 0.25])
set(sp(1),'position',[0.16 0.69 0.8 0.25])
if max(size(ni)) == 1
    leg = legend(pl(1),upper(fn1{1}));
else
    leg = legend(pl,upper(fn1));
end
set(leg,'position',[0.81 0.33 0.05 0.05])
prettyfigures('text',12,'labels',14,...
    'box',1)

%Save Figures!
h = findobj('type','figure');
if max(size(h)) == 2;
    export_fig(figure(h(2)),[sdir expname '_' fn1{1} 'AmpBD'],'-png')
    export_fig(figure(h(1)),[sdir expname '_' fn1{1} 'BDVelSigma_b'],'-png')
else
    export_fig(figure(h(4)),[sdir expname '_' fn1{1} 'AmpBD'],'-png')
    export_fig(figure(h(3)),[sdir expname '_' fn1{2} 'AmpBD'],'-png')
    export_fig(figure(h(2)),[sdir expname '_' fn1{3} 'AmpBD'],'-png')
    export_fig(figure(h(1)),[sdir expname '_BDVelSigma_b'],'-png')
end