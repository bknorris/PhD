%Calculate basic statistics on number of events, wg & ig, acc & ero
clear,close all
load('d:\Projects\Mekong_W2015\DataAnalysis\Paper3\Wavelet\CombEventdata_flood.mat');
data.flood = dat;
load('d:\projects\Mekong_W2015\DataAnalysis\Paper3\Wavelet\CombEventdata_ebb.mat');
data.ebb = dat;clear dat
fn = fieldnames(data);
sfdir = 'f:\GradSchool\DataAnalysis\Paper3\Figures\';
flds = {'mud';'fringe';'forest'};
band = {'wave';'ig'};
nevents = [3 11 2;3 7 3];
nexps = [3 7 2;3 5 3];
for i = 1:2
    disp([fn{i} ' tide'])
    for j = 1:length(flds)
        disp(['Location: ' flds{j}])
        nwave = numel(data.(fn{i}).(flds{j}).wave.eventl);
        nig = numel(data.(fn{i}).(flds{j}).ig.eventl);
        nwvig = nwave/nig;
        fprintf('Number of wave events: %0.1f\n',nwave/nexps(i,j))
        fprintf('Number of IG events: %0.1f\n',nig/nexps(i,j))
%         fprintf('Ratio of wave/IG events: %0.1f\n',(nwave/nexps(i,j))/nig/nexps(i,j))
        aevents = length(find(data.(fn{i}).(flds{j}).wave.deltbd>0))/length(find(data.(fn{i}).(flds{j}).wave.deltbd<0));
        eevents = length(find(data.(fn{i}).(flds{j}).ig.deltbd>0))/length(find(data.(fn{i}).(flds{j}).ig.deltbd<0));
        fprintf('Ratio of accretion/erosion events, wave: %0.1f, and IG: %0.1f\n',aevents,eevents)
        id1 = find(data.(fn{i}).(flds{j}).wave.deltbd>0);id2 = find(data.(fn{i}).(flds{j}).ig.deltbd>0);
        awbdmax = max(data.(fn{i}).(flds{j}).wave.deltbd(id1));aibdmax = max(data.(fn{i}).(flds{j}).ig.deltbd(id2)); %#ok<FNDSB>
        fprintf('Ratio of accretionary maximum BLE ratio, wave/IG events: %0.1f\n',abs(awbdmax/aibdmax))
        id1 = find(data.(fn{i}).(flds{j}).wave.deltbd<0);id2 = find(data.(fn{i}).(flds{j}).ig.deltbd<0);
        ewbdmax = min(data.(fn{i}).(flds{j}).wave.deltbd(id1));eibdmax = min(data.(fn{i}).(flds{j}).ig.deltbd(id2));
        fprintf('Ratio of erosional minimum BLE ratio, wave/IG events: %0.1f\n',abs(ewbdmax/eibdmax))
        %Calculate simple regression of eventl and deltbd for acc & ero
        %events
        eventl = [data.(fn{i}).(flds{j}).wave.eventl; data.(fn{i}).(flds{j}).ig.eventl];
        deltbd = [data.(fn{i}).(flds{j}).wave.deltbd; data.(fn{i}).(flds{j}).ig.deltbd];
        id = find(deltbd>0);
        pf = polyfit(eventl(id),deltbd(id),1);
        yfit = polyval(pf,eventl(id));
        yresid = deltbd(id)-yfit;
        SSresid = sum(yresid.^2);
        SStotal = (length(deltbd(id))-1)*var(deltbd(id));
        rsq = 1-SSresid/SStotal;
        fprintf('R-squared for event length/delta bd, accretion: %0.4f\n',rsq)
        id = find(deltbd<0);
        pf = polyfit(eventl(id),deltbd(id),1);
        yfit = polyval(pf,eventl(id));
        yresid = deltbd(id)-yfit;
        SSresid = sum(yresid.^2);
        SStotal = (length(deltbd(id))-1)*var(deltbd(id));
        rsq = 1-SSresid/SStotal;
        fprintf('R-squared for event length/delta bd, erosion: %0.4f\n',rsq)
    end
end
fprintf('End of basic statistics from event duration file\n\n')
fprintf('Calculating net bed change for all experiments...\n')
dpath = '\DataAnalysis\Paper3\';
ypath = {'Mekong_F2014';'Mekong_W2015'};
%load run file to tell program which files & VPs to load
rdir = 'g:\GradSchool\DataAnalysis\Paper3\ExperimentalDesign\';
filename = {'AutoMLRrunfile_flood.csv';'AutoMLRrunfile_ebb.csv'};
for j = 1:2
    fid = fopen([rdir filename{j}]);
    rfile = textscan(fid,'%s%s%s','delimiter',',');
    rdate = rfile{1};rexp = rfile{2};rfld = rfile{3};
    for i = 1:2
        npath = ['f:\' ypath{i} dpath];
        folders = dir([npath '\BottomTrack\']);
        folders = {folders(3:end).name};
        for ii = 1:length(folders)
            if any(strcmp(folders{ii},rdate))
                vpid = strcmp(folders{ii},rdate);vpid = find(vpid);
                bfile = dir([npath 'BottomTrack\' folders{ii} '\' '*.mat']);
                bd = load([npath 'BottomTrack\' folders{ii} '\',bfile.name]);
            else
                continue
            end
            ffn = fieldnames(bd);
            for k = 1:length(vpid)
                bdist = bd.(ffn{k}).bdist;
                id1 = find(~isnan(bdist),1,'first');
                id2 = find(~isnan(bdist),1,'last');
                dbdist = bdist(id1)-bdist(id2);
                %save to structure
                data.(fn{j}).(rfld{vpid(k)}).dbdist(vpid(k)) = dbdist;
            end
        end
    end
end
for i = 1:2
    disp([fn{i} ' tide'])
    for ii = 1:3
        disp(flds{ii})
        a = data.(fn{i}).(flds{ii}).dbdist;
        a(a==0) = [];
        fprintf('Mean change in BLE over all experiments: %0.4f mm\n',mean(a)*1000)
        fprintf('Median BLE over all experiments: %0.4f mm\n',median(a)*1000)
        fprintf('Max BLE over all experiments: %0.4f mm\n',max(a)*1000)
        fprintf('Min BLE over all experiments: %0.4f mm\n',min(a)*1000)
    end
end
        