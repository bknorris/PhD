%Plotting script for vertical velocity time series and spectra for FSS3 and
%DPS2 downward-looking ADCPs. Spectra will be plotted from velocity data
%(Z).

%in the FSS, AD5116 was deployed above a bare-bed, AD5117 was deployed
%above pneumatophores

clear
close all
ddir = '/Volumes/MEKONG1/Mekong_W2015/Data/Aquadopp/FSS/';
fdir = '/Volumes/MEKONG1/Mekong_W2015/Figures/Spectra&Waves/';

load([ddir 'AD5116_9March2015_f_pad.mat']);
data.a6 = aqdp; clear aqdp
load([ddir 'AD5117_9March2015_f_pad.mat']);
data.a7 = aqdp; clear aqdp
start = datenum(2015,03,07,14,20,00);
stop = datenum(2015,03,07,16,00,00);

%crop data to desired time interval
fn = {'a6';'a7'}; %field names of data
DATA = struct(fn{1},[],fn{2},[]);        
nlin = 1E2;
maxgaps = 1E5;
%copy only the information that's needed from the data files & bridge gaps
%with cmgbridge
for i = 1:2
    ind = find(data.(fn{i}).datenum >= start & data.(fn{i}).datenum <=stop);
    DATA.(fn{i}).dt = data.(fn{i}).datenum(ind);
    DATA.(fn{i}).rb = data.(fn{i}).rangebins;
    DATA.(fn{i}).ncells = data.(fn{i}).nCells;
    DATA.(fn{i}).cellsize = data.(fn{i}).metadata.cellsize;
    DATA.(fn{i}).blankdist = data.(fn{i}).metadata.blankdist;
    DATA.(fn{i}).hab = data.(fn{i}).metadata.HAB/1000; %mm to m
    DATA.(fn{i}).fs = str2double(data.(fn{i}).metadata.samprate(1:2));
    tmp = data.(fn{i}).w(ind,:);
    [m,n] = size(tmp);
    DATA.(fn{i}).vz = zeros(m,n);
    for j = 1:n
        if cmgidgaps(tmp(:,j)) > 0
            disp(['Found ' num2str(cmgidgaps(tmp(:,j))) ' gaps in velocity time-series ' fn{i} ' bin ' num2str(j)])
            tmp2 = cmgbridge(tmp(:,j),nlin,maxgaps,maxgaps);
            if cmgidgaps(tmp2) > 0
                disp(['Number of gaps in velocity remaining: ' num2str(cmgidgaps(tmp2))])
            else
                disp('Gaps filled')
            end
            DATA.(fn{i}).vz(:,j) = tmp2;
        else
            DATA.(fn{i}).vz(:,j) = tmp(:,j);
        end
        clear tmp2
    end
    clear tmp
end
clear data %dump excess


% for i = 1:2
%     %Calculate PSD(timeseries)
%     psvz = detrend(vz);
%     n = length(psvz);    
%     nf = DATA.(fn{i}).fs/2; %nyquist criterium
%     nfft = n/2;
%     [psd,freq] = pwelch(psvz,hanning(n),n/2,n,DATA.(fn{i}).fs);
%     pf = psd(2:end,1);
%     pf = runningmax(pf,100);
%     pf = fastsmooth(pf,50,2,1); %make psd nicer for plotting
%     freq(1) = [];
%     DATA.(fn{i}).psd = pf;
%     DATA.(fn{i}).freq = freq;
% end


