%check the aquadopp directions with the vector from F2F2 and F2F3
clear
aqdpdir = 'c:\Users\bkn5\Projects\Mekong_W2015\Data\Aquadopp\F2F2\';
vecdir = 'c:\Users\bkn5\Projects\Mekong_W2015\Data\Vector\F2F2\';

aqdps = dir([aqdpdir '*.mat']);
%sort dir by date, use only the most recent files
dd = strfind({aqdps.date},'09-Feb-2016'); %date the files were reprocessed
id = find(not(cellfun('isempty', dd)));

%run through all the aquadopps
for i = 1:length(id)
    disp(['Loading ' aqdps(id(i)).name])
    load([aqdpdir aqdps(id(i)).name])
    
    %need to depth then time average, then calculate quivers at a given
    %time interval... say 10 minutes.
    
    lat = aqdp.metadata.lat;
    lon = aqdp.metadata.lon;
    
    u = nanmean(aqdp.u,2);
    v = nanmean(aqdp.v,2);
    
    intv = 10; %minutes
    fs = 8;%hz
    avt = fs*60*intv;
    
    idx = [1 avt:avt:length(u)];
    for j = 1:length(idx)-1
        t = idx(j):idx(j+1)
        U = u(t);
        V = v(t);