%Calculate bottom tracking time series for the F2F2 experiment
%where the VPs are closest to the bed.

%Data files should match the timestep of the companion file; turb.
%dissipation (30 sec intv). 

clear
close all
datdir = 'D:\Projects\Mekong_W2015\Data\Vectrino\7March2015\';
fname = dir([datdir '*070315.mat']);
inst = cell(3,1);
for i = 1:3
    inst{i} = [datdir char({fname(i).name})];
end

ts = 'TKE';
win = 3000/5;
step = 1500/5;

for i = 1:3
    start = datenum(2015,03,07,13,36,00);
    stop = datenum(2015,03,12,17,08,00); %TKE stop time
    vname = char(lower(regexp(inst{i},'(?=\VP)(.*?)(?=\_)','match')));
    
    disp(['Loading ' inst{i}])
    load(inst{i})
    bd = VPRO.Data.BottomCheck_BottomDistance;
    gmt2ict = datenum(0,0,0,1,0,0)*7;
    time = VPRO.Data.BottomCheck_HostTimeMatlab+gmt2ict;
    id = find(time >= start & time <= stop);
    bd = bd(id);time = time(id);
    
    %%%Smooth bottom distance
    bdmid = my_running_median(bd,100);
    %window to smooth for plotting
    idx = [1 step:step:length(time)];
    bdmax = zeros(length(idx),1);
    time2 = zeros(length(idx),1);
    for ii = 1:length(idx)-1
        if abs(length(bdmid)-idx(ii)) < win
            continue
        else
            bwin = bdmid(idx(ii):idx(ii)+win);
            bdmax(ii,:) = max(bwin);
            time2(ii,:) = time(idx(ii));
        end
    end
    bdmax(bdmax == 0) = []; %remove trailing zeros
    time2(time2 == 0) = [];
    %     id2 = find(bdmax < mean(bdmax)-2.5*std(bdmax));
    %     if ~isempty(id2)
    %         if id2(1) < 21
    %             rm = min(id2):max(id2)+20;
    %         elseif max(id2) > length(bdmax)-21
    %             rm = min(id2)-20:max(id2);
    %         else
    %             rm = min(id2)-20:max(id2)+20;
    %         end
    %         bdmax(rm) = NaN;
    %         bad=isnan(bdmax);
    %         gd=find(~bad);
    %         bad([1:(min(gd)-1) (max(gd)+1):end])=0;
    %         bdmax(bad)=interp1(gd,bdmax(gd),find(bad),'pchip');
    %     end
    bdavs = smooth(bdmax,25);
    data.(vname).time = time2;
    data.(vname).bd = bdavs;
    clear VPRO
end
dirname = 'd:\Projects\Mekong_W2015\DataAnalysis\Paper1\';
save([dirname 'BDtrack_' ts],'data','-v7.3')