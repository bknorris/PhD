clear
maindir = '/Volumes/MEKONG3/Mekong_F2014/DataAnalysis/Paper2/';
savedatdir = '/Volumes/MEKONG3/Mekong_F2014/DataAnalysis/Paper2/SNR/';
files = dir([maindir '*Vels.mat']);files = {files.name};
toprocess = 1:8;

for k = 1:length(toprocess)
    disp(['Loading ' files{toprocess(k)}])
    load([maindir files{toprocess(k)}])
    
    fn = fieldnames(dat);
    
    %%%%Averaging settings
    window = 30; %30 second averaging interval
    step = 10; %10 second step
    fs = 50;
    avt = step*fs; %samples/step
    nwin = window*fs; %samples/window
    for i = 1:3                                                                 %loop through instruments VP1-VP3
        disp(['Analysing ' fn{i}])
        nsamp = length(dat.(fn{i}).time);
        ind = [1 avt:avt:nsamp];

        %%%%
        for ii = 1:length(ind)                                                  %loop through time, windowed to the settings
            if abs(nsamp-ind(ii)) < nwin                                        %skip the last few indexes approaching the end of the t-s
                continue
            else
                idx = ind(ii):ind(ii)+nwin-1;
                SNR1 = dat.(fn{i}).SNR1(idx,:);
                SNR2 = dat.(fn{i}).SNR2(idx,:);
                SNR3 = dat.(fn{i}).SNR3(idx,:);
                SNR4 = dat.(fn{i}).SNR4(idx,:);
            end

            %%%% extract time of average
            data.time(ii) = dat.(fn{i}).time(ind(ii));
            
            %%%% average
            for j = 1:35
                data.beam1(ii,j) = mean(SNR1(:,j));
                data.beam2(ii,j) = mean(SNR2(:,j));
                data.beam3(ii,j) = mean(SNR3(:,j));
                data.beam4(ii,j) = mean(SNR4(:,j));
            end
        end
        SNR.(fn{i}) = data;
        SNR.(fn{i}).rb = dat.vpro1.rb;
        SNR.(fn{i}).intv = ['Averaging interval: ' num2str(step) ' seconds'];
    end
    fname = regexprep(files{toprocess(k)},'Vels.mat','');
    ft = [fname 'SNRavgs'];
    save([savedatdir ft],'SNR','-v7.3')
end