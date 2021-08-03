clear
load('d:\Mekong_W2015\DataAnalysis\Paper3\BottomTrack\12-03-15\F2F3_2_bdtrace.mat')
m = matfile('d:\Mekong_W2015\DataAnalysis\Paper3\VPs\12-03-15\12March2015_Vels.mat');
bd = vpro1;clear vpro2 vpro3
dat = m.vpro1;
start = datenum(2015,03,12,7,01,55);
stop = datenum(2015,03,12,9,12,47);

vpt = dat.time;
vpx = dat.x;
vpy = dat.y;
bdt = bd.time;
bds = bd.bdist;
bins = 1:5;
%run pca on x and y to figure out dominant current direction. X points
%north.
out = cmgpca(vpx,vpy,[],0);
heading = mean(out.mdir);
th = heading*pi/180;
R = [cos(th) -sin(th); sin(th) cos(th)];
vpx = nanmean(vpx(:,bins),2);
vpy = nanmean(vpy(:,bins),2);
vpxy = [vpx vpy];rxy = zeros(size(vpxy));
for jj = 1:length(vpxy)
    rxy(jj,:) = vpxy(jj,:)*R;end %rotate x, y to cross & along-shore
%NOTE: some files have long data gaps. Files are cropped to
%user-specified lengths in the 2-3 column of the run file (.csv).
tid = find(vpt>=start&vpt<=stop);
rxy = rxy(tid,:);
vpt = vpt(tid,:);
tid = find(bdt>=start&bdt<=stop);
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
time = (dt:dt:length(x10)*dt)/60; %time in minutes for x-axis plotting
xt = [time' x10];
yt = [time' y10];
bt = [time' bds];
[Rsq,period,scale,coi,sig95,Wxy,t,dt]=wtc(xt,bt,'Dj',(1/10),'S0',(1/60));
plotwtc(Rsq,period,coi,sig95,t,Wxy,dt)

%Plot the results
cb = brewermap(100,'YlOrRd');
colormap(gca,cb)
xlabel('Time (minutes)'),ylabel('Period (minutes)')
cbl = colorbar('peer',gca,'Tag','colorbar1');
ylabel(cbl,'Magnitude Squared Coherence')
