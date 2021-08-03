%Load and reorganize RBR datafiles for easy plotting
clear
dirc = 'C:\Users\bkn5\Projects\Mekong_F2014\Data\RBR\SoloT\Control\';

metadata.deployment_date = [2014 09 27 12 00 00];
metadata.recovery_date = [2014 09 27 17 30 00];
metadata.inst_type = 'RBR Solo T Logger';  % type of instrument and instrument manufacturer
metadata.data_cmt = 'Mekong 2014 Control Study - Mudflat only'; % any comment
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
cd(dirc)
filelist = dir('*.mat');
files = {filelist.name};
for i = 1:length(files)
    load(files{i}); %stuct name: RBR
    disp(['Reading in file ' files{i}]);
    [n,m]=size(RBR.sampletimes);
    dates = zeros(n,m);
    yrdates = zeros(n,m);
    
    %define metadata structure using instrument info
    metadata.instname = RBR.name;
    metadata.serial = str2num(RBR.name(end-5:end));
    metadata.samprate = RBR.sampleperiod;
    metadata.channelnames = RBR.channelnames;
    metadata.chanunits = RBR.channelunits;
    metadata.coefficients = RBR.coefficients;
    metadata.parameters = RBR.parameters;
    
    %make datetime/yearday vectors from sampletimes
    for ii = 1:length(RBR.sampletimes)
        dd = str2num(RBR.sampletimes{ii}(1:2)); %#ok<*ST2NM>
        mm = str2num(RBR.sampletimes{ii}(4:5));
        yyyy = str2num(RBR.sampletimes{ii}(7:10));
        times = datestr(RBR.sampletimes{ii},'HH:MM:SS.FFF');
        hr = str2num(times(1:2));
        mi = str2num(times(4:5));
        sec = str2num(times(7:12));
        dates(ii) = datenum(yyyy,mm,dd,hr,mi,sec)';
        yrdates(ii) = yearday(yyyy,mm,dd+hr/24+mi/(24*60)+sec/(24*60*60));
    end
    
    %crop data to length of deployment
    deploy = datenum(metadata.deployment_date);
    recover = datenum(metadata.recovery_date);
    ind = find(dates >=deploy & dates <= recover);
    fname = RBR.channelnames;
    fname = regexprep(fname,' ','_');
    for k = 1:length(fname)
        thedata.(fname{k}) = RBR.data(ind,k);
    end
    thedata.Datetime = dates(ind);
    thedata.Yearday = yrdates(ind);
    clear RBR
    RBR = thedata;
    RBR.Metadata = metadata;
    sname = ['RBR_' num2str(metadata.serial)];
    save([dirc sname],'RBR')
    disp([sname ' saved'])
end

