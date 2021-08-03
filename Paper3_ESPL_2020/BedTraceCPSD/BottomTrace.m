%Generate Bed Trace for VPs for spectral analysis
%
% Paper 3: Sediment Motion in Mangroves
%
% This is Version 1.0 of this script
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear
close all

%let's start with F2F2_2 as an example
dir1 = 'e:\Mekong_W2015\DataAnalysis\Paper3\VPs\06-03-15\';
files = {'6March2015_Vels.mat';'6March2015_Sen.mat'};
hab = 0.063;
expname = 'F2F2_2';
fn1 = whos('-file',[dir1 files{1}]);fn1 = {fn1.name};
V = matfile([dir1 files{1}]);
S = matfile([dir1 files{2}]);
% load('d:\Projects\Mekong_W2015\Data\Aquadopp\F2F2\AD5116_9March2015.mat')

%load the data
for ii = 1:3
    vp = V.(fn1{ii});
    se = S.(fn1{ii});
    
    bdt = vp.bdtime;
    bdist = vp.bdist;
    amp = (se.Amp1+se.Amp2+se.Amp3+se.Amp4)./4;
    ampt = se.time;
    vph = hab-vp.rb;
    figure(1) %plot amp
    subplot(121)
    imagesc(ampt,vph,amp'),set(gca,'ydir','normal')
    hold on
    plot(bdt,hab-bdist)
    datetick('x','HH:MM','keepticks','keeplimits')
    xlabel(['Time on ' datestr(bdt(1),'dd-mm-yy')])
    ylabel('Height Above Bottom (m)')
    
    sfiledir = 'g:\Mekong_W2015\DataAnalysis\Paper3\BottomTrack\06-03-15\';
    bds = my_running_median(bdist,50);
    bd = my_running_max(bds,1024);
    subplot(122)
    imagesc(ampt,vph,amp'),set(gca,'ydir','normal')
    hold on
    plot(bdt,hab-bd)
    datetick('x','HH:MM','keepticks','keeplimits')
    xlabel(['Time on ' datestr(bdt(1),'dd-mm-yy')])
    pause
    %manually despike timeseries with ginput
    yes = 0;
    while yes == 0
        f1 = figure(2);
        set(f1,'PaperOrientation','portrait',...
            'position',[400 400   1200   400]);
        p(1) = plot(bdt,bd,'k');
        datetickzoom('x','HH:MM:SS','keepticks','keeplimits')
        ylabel('BLE (m)')
        title(['Bottom Trace - ' fn1{ii}])
        hold on
        
        disp('USER select times to despike')
        disp('Use left click to zoom, right click to mark points')
        disp('press RETURN when finished')
        [xp,yp] = ginput_zoom;
        l = length(xp);
        if mod(l,2) ~= 0
            disp('Number of selections must be even')
            close(f1)
            continue
        end
        xs = zeros(l,1);
        for j = 1:length(xp)
            tmp = abs(bdt-xp(j));
            [~,idx] = min(tmp);
            xs(j) = idx;
        end
        pairs = reshape(1:l,2,l/2);[~,t]=size(pairs);
        for j = 1:t
            idx = xs(pairs(1,j)):xs(pairs(2,j));
            bd(idx) = NaN;
        end
        bd = fixgaps(bd);
        bdt(isnan(bd)) = [];bd(isnan(bd)) = [];
        bd = smooth(bd,5);
        p(2) = plot(bdt,bd,'r');
        legend(p,{'No Despike';'Despike'})
        yes = 1;
        pts = xs;
        btrace.(fn1{ii}).time = bdt;
        btrace.(fn1{ii}).bdist = bd;
        pause(2),close(f1)
    end
end
save([sfiledir expname '_bdtrace_v2'],'-struct','btrace')