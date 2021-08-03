function aqdp = padbursts(a)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Pad Aquadopp burst data with NaNs to account for the seconds in between
%  bursts for burst sampling instruments. This function will also
%  pad the Aquadopp correlations, backscatter, and time variables.
%
%  Usage:
%  aqdp = padbursts(a)
%
%  Where: a is an aquadopp data file, aqdp is the resultant datafile
%
%  Script developed by Benjamin K Norris, University of Waikato, NZ 2015
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if strcmp(a.metadata.samprate,'1 Hz') == 1
    sampmode = 'CONTINUOUS';
    msg = 'The instrument was set to record continuously. Exiting...';
    warning(msg)
    aqdp = a;
else
    sampmode = 'BURST';
    samprate = regexprep(a.metadata.samprate,'Hz','');
    srate = str2num(samprate);
    spb = a.metadata.spb;
    blength = spb/srate; %burst length
    speriod = a.metadata.dt; %sample period
    missec = speriod-blength; %missing seconds
    
    %find the indices of the end of the bursts
    limit = datenum(0,0,0,0,0,missec);
    int = datenum(0,0,0,0,0,missec/(missec*srate)); %convert miliseconds to datenumbers
    diffs = diff(a.datenum);ind = find(diffs > limit);
    snames = fieldnames(a);
    dat = struct();
    idx = [0; ind; numel(a.datenum)];
    
    %Pad date variables
    for k = 1:length(idx)-1
        daten = a.datenum(idx(k)+1:idx(k+1));
        yeard = a.yearday(idx(k)+1:idx(k+1));
        f1 = idx(k+1);f2 = length(daten);
        for j = 1:(missec*srate)
            dfill(j,:) = a.datenum(f1)+(int*j);
            yfill(j,:) = a.yearday(f1)+(int*j);
            daten(f2+j,:) = dfill(j,:);
            yeard(f2+j,:) = yfill(j,:);
        end
        dat.datenum{k,1} = daten;
        dat.yearday{k,1} = yeard;
        
    end
    dat.datenum = cat(1,dat.datenum{:});
    dat.yearday = cat(1,dat.yearday{:});
    %Pad all other variables with NaNs
    
    for ii = [3:19 22 26:34]; %clip all fields to start/stop times
        for k = 1:length(idx)-1
            field = a.(snames{ii})(idx(k)+1:idx(k+1),:);
            [~,w] = size(a.(snames{ii}));
            fill = NaN(missec*srate,w);
            field = [field; fill];
            dat.(snames{ii}){k,1} = field;
        end
        a = rmfield(a,(snames{ii}));
        dat.(snames{ii}) = cat(1,dat.(snames{ii}){:});
        disp([snames{ii} ' bursts padded with NaNs'])
    end
    a = rmfield(a,{'datenum','yearday'});
    aqdp = catstruct(a,dat);
end
end




