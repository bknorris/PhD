%Calculate pressure spectra from an Aquadopp
clear

disp('Calculating Pressure Spectra')
load ab_nu5percent %stored in C:\Users\Documents\MATLAB\
load HR3_12March2015_f_pad.mat

name = 'ADHR3';
savefigdir = 'c:\Users\bkn5\Projects\Mekong_W2015\Figures\Aquadopps\'; 


%crop times to just the day of interest
% figure
% plot(aqdp.datenum,aqdp.pressure)
% datetickzoom('x','dd HH:MM:SS','keepticks','keeplimits')
% pause
startt = datenum(2015,03,10,14,15,00);
endt = datenum(2015,03,10,17,00,00);
ind = find(aqdp.datenum >= startt & aqdp.datenum <= endt);
dat = aqdp.pressure(ind); %#ok<FNDSB>

%averaging interval
intv = 10; %averaging window in minutes (must be >= 1 min)
sr = 8; %sample rate (Hz)
avt = 60*intv*sr;
%define spectra parameters
nf = sr/2; %nyquist criterium
nfft = 0.15*avt;
overlap = 0.7;
overlpts = floor(overlap*nfft);

%need to just find the points that aren't zero
ind = find(dat>0);dat = dat(ind);
dat = detrend(dat);
[n,~] = size(dat);
idx = avt:avt:n;
idxx = [1 idx]; %do not include remainders here for this calculation

%preallocate variables
avdat = zeros(length(idx),1);
%calculate specta of field
k = 1;
while k < length(idxx)
    idx = idxx(k):idxx(k+1);
    for j = 1:length(idx)
        avdat(j,:) = dat(idx(j),:);
    end
    [spec(:,k),f(:,k)] = cpsd(avdat,avdat,hanning(nfft,'periodic'),overlpts,nfft,sr); %#ok<*SAGROW>
    k = k+1;
end
%calculate degrees of freedom and the confidence intervals
dof = floor(8*length(avdat)/(3*nfft));
fprintf('Degrees of freedom: %d\n',dof)
ciLimit = [0.025,0.975]; % For 95% CIs.
ab = ab_nu5(dof,:);
ci_up=dof.*(spec)./ab(1); %Upper
ci_lw=dof.*(spec)./ab(2); %Lower

f1 = figure;
set(f1,'PaperOrientation','portrait',...
    'position',[400 200   1000   800]);
set(gcf,'color','w','PaperPositionMode','auto')
[~,m] = size(spec);
c = jet(m);

for ii = 1:m
    area(f(:,ii),ci_up(:,ii),'FaceColor',[0.9 0.9 0.9],'LineStyle','none'),hold on
    area(f(:,ii),ci_lw(:,ii),'FaceColor',[1 1 1],'LineStyle','none'),hold on
    p(ii) = plot(f(:,ii),spec(:,ii),'-x','linewidth',1.5);
    set(p(ii),'Color',c(ii,:))
    hold on
    box on
    xlabel('\bf\itHz')
    ylabel('\bf\itmbar^2/Hz')
    names{ii} = sprintf('%d min to %d min',(intv*ii-intv),intv*ii);
end
grid on
leg = legend(p,names,'location','northeastoutside');
set(gca,'Xlim',[0.05 1],'Ylim',[0 500],...
            'GridLineStyle',':')
title(['\bf\it' name ' Pressure Spectra, ' datestr(startt,'dd/mm/yyyy')])
% prompt = 'Save Figure? [y/n] ';
% result = input(prompt,'s');
% if strcmp(result,'y') || strcmp(result,'yes');
%     fpath = savefigdir;fname = [name '_PresSpectra'];
%     export_fig([fpath fname],'-png','-m1','-r900','-opengl')
%     disp(['Figure ' fname '.png saved'])
% end
% if strcmp(result,'n') || strcmp(result,'no');
% end

