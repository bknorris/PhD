clear
%define working paths
dpath = '\DataAnalysis\Paper3\';
ypath = {'Mekong_F2014';'Mekong_W2015'};
%load run file to tell program which files & VPs to load
rdir = 'f:\MemoryStick\GradSchool\DataAnalysis\Paper3\ExperimentalDesign\';
fid = fopen([rdir 'AutoTKErunfile_v4.csv']);
rfile = textscan(fid,'%s%s%s%s%n%n%n','delimiter',',');
rexpt = rfile{1};rdate = rfile{2}; %experiment name, experiment date
rstart = rfile{3};rstop = rfile{4}; %start time, stop time (based on VPs)
rinst = rfile{5};rhab = rfile{6}; %instrument name in folder, instrument height above bed
rd = rfile{7}; %mean stem diameter
start = datenum(strcat(rdate,{' '},rstart),'dd-mm-yy HH:MM:SS');
stop = datenum(strcat(rdate,{' '},rstop),'dd-mm-yy HH:MM:SS');
stop = stop+datenum(0,0,0,0,3,0); %account for windowed indexing
for i = 1:2
    npath = ['f:\' ypath{i} dpath 'BottomTrack\'];
    folders = dir(npath);
    folders = {folders(3:end).name};
    for ii = 1:length(folders)
        if any(strcmp(folders{ii},rdate))
            inid = strcmp(folders{ii},rdate);
            inid = find(inid);
            tpath = ['f:\' ypath{i} dpath 'Turbulence\' folders{ii} '\'];
            bpath = [npath folders{ii} '\'];
            disp(['Loading files from ' tpath])
            disp(['Loading files from ' bpath])
        else
            continue
        end
        tfile = dir([tpath,'*_TKE.mat']);
        bfile = dir([bpath,'*bdtrace_v2.mat']);
        tke = load([tpath tfile.name]);
        bd = load([bpath bfile.name]);
        fn = fieldnames(tke);
        f(ii) = figure(ii);
        c = [207 176 126;
            60 166 74;
            4 76 41]./255;
        for j = 1:length(fn)
            bh = bd.(fn{j}).bdist(1)-bd.(fn{j}).bdist;
            bdt = bd.(fn{j}).time;
            tdt = tke.(fn{j}).time;
            eps = (tke.(fn{j}).z1.E+tke.(fn{j}).z2.E)./2;
            [~,idx] = unique(tdt);
            ep10 = zeros(30,length(bdt));
            for k = 1:30
                ep10(k,:) = interp1(tdt(idx),eps(k,idx),bdt,'linear','extrap');
            end
            %04/17/20: Find depth bins that are a consistent elevation
            %above bed level for TKE extraction!
            vh = bd.(fn{j}).bdist(1)-0.04-linspace(0.001,0.03,30);
            eps = zeros(length(bdt),1);
            zed = zeros(length(bdt),1);
            for k = 1:length(bdt)
                [minval,bid] = min(abs((bh(k)+0.007)-vh));
                eps(k) = ep10(bid,k);
                zed(k) = vh(bid);
            end
            % Convert Epsilon to Turbulence intensity (k)
            d = rd(inid(j));
            if d == 0 %this is a mudflat experiment
                C0 = 0.19;
                kappa = 0.41;
                TKE = ((eps.^(2/3)).*((kappa*abs(zed)).^(2/3)))./C0;
            elseif d > 0 %this is a vegetated experiment
                alpha = 0.9; 
                TKE = 2*(eps.^(2/3))*d^(2/3)*alpha^2;
            end

            %Plot routine
            sp(1) = subplot(211);
            if j == 1
                plot(bdt,zeros(size(bdt)),':k','linewidth',1.5),hold on
            end
            plot(bdt,(bh).*1000,'color',c(j,:),'linewidth',1.5)
            sp(2) = subplot(212);
            p(j) = plot(bdt,eps,'color',c(j,:),'linewidth',1.5); hold on
        end
        set([sp(1) sp(2)],'xlim',[min(bdt) max(bdt)])
        set(sp(1),'xticklabel',[])
        datetick('x','HH:MM','keepticks','keeplimits')
        title(sp(1),[rdate{inid(j)}])
        ylabel(sp(1),'Height Above Bottom [mm]')
        ylabel(sp(2),'TKE Dissipation Rate [W/m^2]')
%         ylabel(sp(2),'k [m^2/s^2]')
        xlabel(sp(2),'Time of Day [HH:MM]')
        leg = legend(p,{'Mud';'Fringe';'Forest'},'position',[0.765 0.4 0.15 0.15]);
        
        prettyfigures('text',10,'labels',12,'box',1)
%         sfpath = 'e:\Mekong_W2015\DataAnalysis\Paper3\Turbulence\Figures\';
        sfpath = 'e:\MemoryStick\GradSchool\DataAnalysis\Paper3\Figures\';
%         export_fig([sfpath rdate{inid(1)} '_bedlevelTKEts'],'-pdf')
    end
end