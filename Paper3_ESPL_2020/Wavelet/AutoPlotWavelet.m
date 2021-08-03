%This program automatically loads and plots the cross-spectral wavelet analysis.
clear
%define working paths
dpath = '\DataAnalysis\Paper3\';
ypath = {'Mekong_F2014';'Mekong_W2015'};
%load run file to tell program which files & VPs to load
rdir = 'e:\MemoryStick\GradSchool\DataAnalysis\Paper3\ExperimentalDesign\';
fid = fopen([rdir 'AutoWLrunfile_v6.csv']);
rfile = textscan(fid,'%s%s%s%n%n%n','delimiter',',');
rdate = rfile{1};
cb = brewermap(100,'YlOrRd');
for i = 2
    npath = ['e:\' ypath{i} dpath];
    folders = dir([npath '\Wavelet\']);
    folders = {folders(3:end).name};
    for ii = 1:length(folders)
        if any(strmatch(folders{ii},rdate))
            vpid = strmatch(folders{ii},rdate);
            disp(['Loading files from ' npath 'Wavelet\' folders{ii} '\'])
            vfile = dir([npath 'Wavelet\' folders{ii} '\','*_wvlt.mat']);
            vfile = {vfile.name};
        else
            continue
        end
        for j = 1:length(vfile)
            load([npath 'Wavelet\' folders{ii} '\' vfile{j}])
            vname = regexprep(vfile{j},'_wvlt.mat','');
            %plot wavelet x-bd then y-bd
            figure(1)
            plotwtc(wvlt.x.Rsq,wvlt.x.period,wvlt.x.coi,wvlt.x.sig95,wvlt.x.t,wvlt.x.Wxy,wvlt.x.dt)
            colormap(gca,cb)
            xlabel('Time (minutes)'),ylabel('Period (minutes)')
            cbl = colorbar('peer',gca,'Tag','colorbar1');
            ylabel(cbl,'Magnitude Squared Coherence')
            title([upper(vname) ' X-BDT Wavelet Coherence ' folders{ii}])
            %%%
            figure(2)
            plotwtc(wvlt.y.Rsq,wvlt.y.period,wvlt.y.coi,wvlt.y.sig95,wvlt.y.t,wvlt.y.Wxy,wvlt.y.dt)
            colormap(gca,cb)
            xlabel('Time (minutes)'),ylabel('Period (minutes)')
            cbl = colorbar('peer',gca,'Tag','colorbar1');
            ylabel(cbl,'Magnitude Squared Coherence')
            title([upper(vname) ' Y-BDT Wavelet Coherence ' folders{ii}])
            %save figs
            figs = findall(allchild(0),'type','figure');
            set(figs(1),'renderer','opengl')
            spath = [npath 'Wavelet\' folders{ii} '\'];
            prettyfigures('text',13,'labels',14,'box',1) 
            export_fig(figs(1),[spath vname '_ybd'],'-png')
            export_fig(figs(2),[spath vname '_xbd'],'-png')
            close all
            clear wvlt
        end
    end
end
            