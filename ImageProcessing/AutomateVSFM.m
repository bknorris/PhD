%run the whole visualSFM suite (pair matching & pmvs) on the 'ToRun'
%quadrats

clear
path1 = 'C:\Users\bkn5\Projects\Mekong_W2015\Images\Quadrats\ToRun\';
dirs = regexp(genpath(path1),['[^;]*'],'match');
tic
for i = 2:length(dirs)
    cd(dirs{i})
    %do the visualsfm stuff, according to Jean...
    filep = [dirs{i} '\'];
    file = 'DenseRec.nvm';
    nn = [filep file];
    [status,cmdout] = system(['C:\VisualSFM\VisualSFM sfm+shared+pairs+pmvs ' filep ' ' nn],'-echo'); 
    %echo should print output to matlab; format here must be
    
    if status == 0
        fprintf('\n')
        disp(['Run ' num2str(i) ' completed; ' dirs{i}])
    else
        warning(['Run ' num2str(i) ' NOT completed; ' dirs{i}])
    end
    
end
fprintf('\nALL 2015 IMAGE QUADRATS COMPLETE \n')
disp([num2str(toc/60) ' minutes elapsed'])
disp('Beginning 2014 image processing now...')
disp(['Run started: ' datestr(now)])
path2 = 'c:\Users\bkn5\Projects\Mekong_F2014\Images\Quadrats\ToRun\';
dirs = regexp(genpath(path2),['[^;]*'],'match');
rstart = tic;
for i = 2:length(dirs)
    start = tic;
    cd(dirs{i})
    %do the visualsfm stuff, according to Jean...
    filep = [dirs{i} '\'];
    file = 'DenseRec.nvm';
    nn = [filep file];
    [status,cmdout] = system(['C:\VisualSFM\VisualSFM sfm+shared+pairs+pmvs ' filep ' ' nn],'-echo'); 
   
    if status == 0
        disp(['Run ' num2str(i) ' completed; ' dirs{i}])
    else
        warning(['Run ' num2str(i) ' NOT completed; ' dirs{i}])
    end
    disp(['Run ' num2str(i) ' completed in ' num2str(toc(start)/60) ' minutes'])
end
fprintf('\nALL 2014 IMAGE QUADRATS COMPLETE \n')
disp(['Total time elapsed: ' num2str(toc(rstart)/60) ' minutes'])
fprintf('\nIMAGE PROCESSING COMPLETE \n')