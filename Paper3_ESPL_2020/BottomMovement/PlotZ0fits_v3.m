clear, close all
load('E:\Mekong_W2015\DataAnalysis\Paper3\VPs\06-03-15\6March2015_Vels.mat')
clear vpro2 vpro3
%load velocities, rotate to heading
rdate = '06-03-15';rstart = '13:37:46';rstop = '17:01:40';
start = datenum(strcat(rdate,{' '},rstart),'dd-mm-yy HH:MM:SS');
stop = datenum(strcat(rdate,{' '},rstop),'dd-mm-yy HH:MM:SS');
stop = stop+datenum(0,0,0,0,3,0); %account for windowed indexing
%crop data to start/stop times
vid = find(vpro1.time>=start&vpro1.time<=stop);
vpx = vpro1.x(vid,:);vpy = vpro1.y(vid,:);
vph = 0.0667; %height of VecPro from bottom track data
heading = -20;
th = heading*pi/180;
R = [cos(th) -sin(th); sin(th) cos(th)];
xr = zeros(size(vpx));yr = zeros(size(vpy));
for k = 1:length(vpx)
    xx = vpx(k,:);
    yy = vpy(k,:);
    for kk = 1:length(xx)
        rxy = [xx(kk) yy(kk)]*R;
        xr(k,kk) = rxy(1);
        yr(k,kk) = rxy(2);
    end
end
%Set up averaging window (we will just take one interval for this example)
win = 180; %seconds
step = 1; %seconds
fs = 50;
avt = fs*step;
nwin = fs*win;
ind = [1 avt:avt:length(vpro1.time)];
time = zeros(length(ind),1);
u = zeros(length(ind),35);zuv = vph-vpro1.rb;c = 1;
for jj = 1:length(ind)
    if abs(length(xr)-ind(jj)) < nwin  %skip the last few indexes approaching the end of the t-s
        continue
    else
        idx = ind(jj):ind(jj)+nwin-1;
    end
    u(c,:) = nanmean(sqrt(xr(idx,:).^2));
    time(c) = vpro1.time(idx(1));
    c = c+1;
end
tid = 1203; %timestep of interest
disp(datestr(time(tid)))
f1 = figure(1);
set(f1,'PaperOrientation','portrait',...
    'position',[400 100   500   400],...
    'renderer','painters');hold on
umean = u(tid,1:25); %cut off bottom
z = zuv(1:25);
%line of best fit
pf = polyfit(umean(19:23),log(z(19:23)),1);
xx = linspace(0.005,0.07,10);
pv = polyval(pf,xx);
R = corrcoef(umean(19:23),log(z(19:23)));
Rsq = R(1,2).^2;fprintf('r-squared: %0.2f\n',Rsq)
%highlight points of BFL
plot(umean,log(z),'.',...
    'markersize',7,...
    'color','k'),hold on
plot(xx,pv,'-r')

title('Site 1, March 6^{th}, 2015')
ylabel('ln(z) [m]')
xlabel('u [m/s]')
prettyfigures('text',12,'labels',13,'box',1,'grid',1,'gstyle',':')
sfdir = 'c:\Users\Bnorr\Documents\GradSchool\Writing\GM_2018_Norris\Figures\V6\';
% export_fig([sfdir 'z0_fit_v3'],'-pdf','-nocrop')
