%Plot spectral slope extraction method using VP1 of HTA1.
clear
load('D:\Projects\Mekong_W2015\DataAnalysis\Paper3\VPs\7March2015_Vels.mat')
heading = 20;
bins = 1:5;
ind = 1:length(dat.vpro1.time)-1;
time = dat.vpro1.time(ind);
x = dat.vpro1.x(ind,bins);
y = dat.vpro1.y(ind,bins);
z = (dat.vpro1.z1(ind,bins)+dat.vpro1.z1(ind,bins))./2;
%I am interested in seeing what the slope of velocity differences look like
%Calculate difference between a bin in the middle of the profile and an
%adjacent bin
mid = median(bins);
xm1 = x(:,mid+1);
xm2 = x(:,mid);
xdif = xm1-xm2;
x = nanmean(x,2);y = nanmean(y,2);z = nanmean(z,2);
rot = (pi*heading)/180;
T = [cos(rot) -sin(rot);...
    sin(rot) cos(rot)];
vels = [x y];
V = vels*T';
xr = V(:,1);yr = V(:,2);
clear dat
%loop in time
fs = 50;
avt = fs*60;
nwin = fs*300;
swin = fs*20; %10 second averaging window (to smooth)
DOF = round((nwin/swin)*2);
nsamp = length(time);
ind = [1 avt:avt:nsamp];

%Pick a moment in time
datdir = 'd:\Projects\Mekong_W2015\DataAnalysis\Paper2\OldFiles\';
load([datdir 'SpecSlopes (2).mat'])
Us = Spec.day1.vpro1.Uslope;Ts = Spec.day1.vpro1.time;
% figure
% plot(Ts,Us,'k'),datetickzoom('x','HH:MM:SS','keepticks','keeplimits')
% [~,id] = min(abs(Ts-b1.Position(1))); %selected point
id = 84;

idx = ind(id):ind(id)+nwin-1;

u = xr(idx);
v = yr(idx);
w = z(idx);
xdif = xdif(idx);
time2 = time(ind(id));

[Cuu,F] = pwelch(u,hanning(swin),swin*0.5,nwin,fs); %cross shore
[Cvv,~] = pwelch(v,hanning(swin),swin*0.5,nwin,fs); %along shore
[Cww,~] = pwelch(w,hanning(swin),swin*0.5,nwin,fs); %vertical
[Cxd,~] = pwelch(xdif,hanning(swin),swin*0.5,nwin,fs); %x-diff
%             figure(1)
%             loglog(F,Cvv,'color',c(j,:)),hold on

minf = 8.2;
maxf = 17;
fc = find(F >= minf & F <= maxf);
%cross shore
inert = Cuu(fc);ff = F(fc);
pu = polyfit(log(ff),log(inert),1);
Uslope = pu(1);
%along shore
inert = Cvv(fc);
pv = polyfit(log(ff),log(inert),1);
Vslope = pv(1);
%vertical
inert = Cww(fc);
pw = polyfit(log(ff),log(inert),1);
Wslope = pw(1);
fid = find(F >= 0.1 & F<= 25);

%Plotting routine
f1 = figure(1);
set(f1,'PaperOrientation','portrait',...
    'position',[400 100   500   400]);
set(gcf,'color','w','paperpositionmode','auto')   
symb = {'p';'d';'^'};
%U
loglog(F(fid),Cuu(fid),'-','Color',[0.6 0.6 0.6]),hold on
%V
loglog(F(fid),Cvv(fid),'-','Color',[0.5 0.5 0.5])
y_hat=exp(pu(1)*log(ff)+pu(2));
loglog(ff,y_hat,'-k','linewidth',1.5);
y_hat=exp(pv(1)*log(ff)+pv(2));
loglog(ff,y_hat,'-k','linewidth',1.5);
%W
y_hat=exp(pw(1)*log(ff)+pw(2));
loglog(F(fid),Cww(fid),'-','Color',[0.4 0.4 0.4])
loglog(ff,y_hat,'-k','linewidth',1.5);
%Plot reference line
xs = linspace(5,25,length(Cuu(fid)));
ys = 0.01.*(xs.^(-5/3));
plot(xs,ys,'-.','Color','k','LineWidth',1);
text(11,2.5E-4,'^-^5^/^3','FontSize',12)
text(0.06,0.006,'\bf\itS_u_u','FontSize',12)
text(0.06,0.002,'\bf\itS_v_v','FontSize',12)
text(0.06,0.00045,'\bf\itS_w_w','FontSize',12)
hold off
%Plot Adjustments
set(gca,...
    'xlim',[0.05 30],...
    'ylim',[10^-7 10^-1],...
    'Linewidth',1.5)
ylabel('Spectral Density (m^2s^-1)','FontSize',13)
xlabel('f (Hz)','FontSize',13)
title(['HTA1 5-min spectra, ' datestr(Ts(id),'HH:MM')],'FontSize',14)
figdir = 'c:\Users\Bnorr\Documents\GradSchool\Writing\JGR_2018_Norris\Figures\Draft2_Figures\Versions\V2\';
% export_fig([figdir 'SlopeExtMethod'],'-pdf')

%Figure 2: Slope of velocity differences
f2 = figure(2);
set(f2,'PaperOrientation','portrait',...
    'position',[400 100   500   400]);
set(gcf,'color','w','paperpositionmode','auto')   
loglog(F(fid),Cxd(fid),'-','Color',[0.4 0.4 0.4]),hold on
loglog(F(fid),Cuu(fid),'-','Color',[0.6 0.6 0.6]),hold on
xs = linspace(5,25,length(Cuu(fid)));
ys = 0.01.*(xs.^(-5/3));
plot(xs,ys,'-.','Color','k','LineWidth',1);
text(11,2.5E-4,'^-^5^/^3','FontSize',12)
