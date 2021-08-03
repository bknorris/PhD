%WaveletSpectra: Compare the local spectrum from the wavelet time series to
%a cpsd of the two wavelet variables, x-vel and bottom trace. Follow the
%same steps for preprocessing as in AutomatedWavelet.m
clear

%define working paths
dpath = '\DataAnalysis\Paper3\';
ypath = {'Mekong_F2014';'Mekong_W2015'};
%load run file to tell program which files & VPs to load
rdir = 'g:\Mekong_W2015\DataAnalysis\Paper3\';
fid = fopen([rdir 'AutoWLrunfile.csv']);
rfile = textscan(fid,'%s%s%s%n%n%n','delimiter',',');
rdate = rfile{1};rstart = rfile{2};%dates corresponding to folder names; start time for crop
rstop = rfile{3};rvp = rfile{4};  %stop time for crop; vp #
rbin5 = rfile{5};rhead = rfile{6}; %use bins 1-5 for averaging; VP heading for cross & along-shore rot
start = datenum(strcat(rdate,{' '},rstart),'dd-mm-yy HH:MM:SS');
stop = datenum(strcat(rdate,{' '},rstop),'dd-mm-yy HH:MM:SS');
dfn = {'vpro1';'vpro2';'vpro3'};
npath = ['g:\' ypath{2} dpath];
folders = dir([npath '\VPs\']);
folders = {folders(3:end).name};
%%%
vpid = strmatch(folders{1},rdate);
disp(['Loading files from ' npath 'VPs\' folders{1} '\'])
disp(['Loading files from ' npath 'BottomTrack\' folders{1} '\'])
vfile = dir([npath 'VPs\' folders{1} '\','*_Vels.mat']);
bfile = dir([npath 'BottomTrack\' folders{1} '\','*_bdtrace.mat']);
bd = load([npath 'BottomTrack\' folders{1} '\' bfile.name]);
wfile = dir([npath 'Wavelet\' folders{1} '\','*_wvlt.mat']);

%just load the variables you need, saves memory

disp(['Processing ' dfn{rvp(vpid(2))}])
dat = load([npath 'VPs\' folders{1} '\' vfile.name],dfn{rvp(vpid(2))});
load([npath 'Wavelet\' folders{1} '\' wfile(2).name]);
vpt = dat.(dfn{rvp(vpid(2))}).time;
vpx = dat.(dfn{rvp(vpid(2))}).x;
vpy = dat.(dfn{rvp(vpid(2))}).y;
bdt = bd.(dfn{rvp(vpid(2))}).time;
bds = bd.(dfn{rvp(vpid(2))}).bdist;
if rbin5(vpid(2)) == 1 %if VP is close to ground, use bins 1-5
    bins = 1:5;
else
    bins = 13:17;
end
heading = rhead(vpid(2));th = heading*pi/180;
R = [cos(th) -sin(th); sin(th) cos(th)];
vpx = nanmean(vpx(:,bins),2);
vpy = nanmean(vpy(:,bins),2);
vpxy = [vpx vpy];rxy = zeros(size(vpxy));
for jj = 1:length(vpxy)
    rxy(jj,:) = vpxy(jj,:)*R;end %rotate x, y to cross & along-shore
%NOTE: some files have long data gaps. Files are cropped to
%user-specified lengths in the 2-3 column of the run file (.csv).
tid = find(vpt>=start(vpid(2))&vpt<=stop(vpid(2)));
rxy = rxy(tid,:);
vpt = vpt(tid,:);
tid = find(bdt>=start(vpid(2))&bdt<=stop(vpid(2)));
bds = bds(tid);
bdt = bdt(tid);
%Downsample velocity to 10Hz
[vpt,idx]=unique(vpt);
x10 = interp1(vpt,rxy(idx,1),bdt);nid = find(isnan(x10));
disp(['Found ' num2str(length(nid)) ' NaNs in x10'])
x10(nid) = 0;
y10 = interp1(vpt,rxy(idx,2),bdt);nid = find(isnan(y10));
disp(['Found ' num2str(length(nid)) ' NaNs in y10'])
y10(nid) = 0;
nid = find(isnan(bds));
disp(['Found ' num2str(length(nid)) ' NaNs in bottom trace'])
bds(nid) = 0;
%Karin recommends detrending each t-s, Grinstead et al. 2004
%recommends plotting the histograms of each t-s to see if they
%are roughly gaussian. NOTE: I have done this, detrending seems
%to make little difference for the velocities, but a big
%difference (in terms of being closer to gaussian) for the
%bottom trace t-s.
x10 = detrend(x10);y10 = detrend(y10);
bds = detrend(bds);
fs = 10;dt = 1/fs;
%plot the wavelet, extract a time
plotwtc(wvlt.x.Rsq,wvlt.x.period,wvlt.x.coi,wvlt.x.sig95,wvlt.x.t,wvlt.x.Wxy,wvlt.x.dt)
%spectra run from 5 min to 12 min
start = bdt(1)+datenum(0,0,0,0,5,0);
stop = bdt(1)+datenum(0,0,0,0,12,0);
id = find(bdt>=start&bdt<=stop);
xx = x10(id);bb = bds(id);
nwin = round(length(bb)/4);
swin = nwin/2;
[pxx,f]=mscohere(xx,bb,hanning(swin),round(swin*0.7),nwin,fs);
%now extract the local spectrum from the wavelet t-s
id = find(wvlt.x.t>=5&wvlt.x.t<=12);
wxy = real(mean(wvlt.x.sig95(:,id),2))-0.14;
period = wvlt.x.period; %10hz

figure
semilogx(1./f,pxx,'r'),hold on
semilogx(period,wxy)
%lets try calculating a moving spectra
avt = 10;
nsamp = length(x10);
ind = [1 avt:avt:nsamp];
msc = zeros(length(ind),swin+1);
F = zeros(length(ind),swin+1);
tt = zeros(length(ind),1);
for ii = 1:length(ind)                                                  %loop through time, windowed to the settings
    if abs(nsamp-ind(ii)) < nwin                                        %skip the last few indexes approaching the end of the t-s
        continue
    else
        id = ind(ii):ind(ii)+nwin-1;
        xx = x10(id);
        bb = bds(id);
        [pxx,f]=mscohere(xx,bb,hanning(swin),round(swin*0.7),nwin,fs);
        msc(ii,:) = pxx;
        F(ii,:) = f;
        tt(ii) = bdt(id(1));
    end
end
%cannot seem to recreate wavelet with msc. 