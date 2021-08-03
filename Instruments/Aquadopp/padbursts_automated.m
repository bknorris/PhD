%padbursts
clear
dirc = 'c:\Users\bkn5\Projects\Mekong_W2015\Data\Aquadopp\';
folders = dir(dirc);dirs = {folders.name};
dirn = [dirc dirs{5} '\'];

files = dir([dirn '*_f.mat']);
for i = 1:length(files)
    filelist = {files.name};varname = filelist{i};
    if ~exist([regexprep(varname,'.mat','') '_pad.mat'])
        load([dirn varname])
        aqdp = padbursts(aqdp);
        fname = [aqdp.metadata.name '_pad.mat'];
        save([dirn fname],'aqdp')
        disp(['file ' fname ' saved'])
    end
end
disp(dirs)