clear, close all
fdir = 'd:\Mekong_W2015\Data\Aquadopp\F2F2\';
fileList = {'HR3_7March2015.mat';'AD5116_9March2015.mat';'AD5117_9March2015.mat'};
load('D:\Mekong_W2015\DataAnalysis\DataReports\Vegetation\VegDat_by_height_1m.mat')
start = datenum(2015,03,05,15,32,10);
stop = datenum(2015,03,05,16,52,30);
data = struct;
fields = {'mud';'fringe';'forest'};
for i = 1:3
    load([fdir fileList{i}])
    time1 = aqdp.datenum;
    vid = find(time1 >= start & time1 <= stop);
    hp = aqdp.metadata.HAB/1000; %height of pressure sensor
    
    U = nanmean(aqdp.u(vid),2);
    V = nanmean(aqdp.v(vid),2);
    p = nanmean(aqdp.pressure(vid),2)+hp;
    x = (sin(aqdp.metadata.lat/57.29578))^2;
    fs = str2double(regexprep(aqdp.metadata.samprate,' Hz',''));
    g = zeros(length(p),1);h = zeros(length(p),1);
    for j = 1:length(p)
        g(j,:) = 9.780318*(1.0+(5.2788E-3+2.36E-5*x)*x)+1.092E-6*p(j,:);
        h(j,:) = ((((-1.82E-15*p(j,:)+2.279E-10)*p(j,:)-2.2512E-5)*p(j,:)+9.72659)*p(j,:))/g(j,:);
    end
    
    data.(fields{i}).hmean = nanmean(h);
    Umean = runningmean(U,fs*60*20);
    Vmean = runningmean(V,fs*60*20);
    data.(fields{i}).ustd = nanstd(U-Umean);
    data.(fields{i}).vstd = nanstd(V-Vmean);
    if i == 2
        data.(fields{i}).a = veg.five.Q1_1A.a;
        data.(fields{i}).z = veg.five.Q1_1A.z';
    elseif i == 3
        data.(fields{i}).a = veg.five.Q1_4A.a;
        data.(fields{i}).z = veg.five.Q1_4A.z';
    else
        data.(fields{i}).a = [];
        data.(fields{i}).z = [];
    end
end
tdir = 'd:\Mekong_F2014\Data\SergioElevationPoints\';
elev = csvread([tdir 'SWTransect.csv']);
data.transect.dist = elev(:,1);
data.transect.elev = elev(:,2);
save(['d:\Mekong_W2015\DataAnalysis\Paper3\' 'DataForSteve'],'data')

