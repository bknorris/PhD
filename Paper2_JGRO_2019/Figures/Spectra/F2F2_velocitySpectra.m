%Generate vertical velocity Spectra from F2F2 for each velocimeter along the
%transect, ADVs, ADCPs, and Vectrinos. As some of these are collocated,
%decide which to plot using a switch.
clear
close all

%%%Time Settings%%%
% start = datenum(2015,03,06,13,50,00);   %low
% stop = datenum(2015,03,06,14,00,00);
% titl = 'Low Tide';
% start = datenum(2015,03,06,14,35,00);     %mid
% stop = datenum(2015,03,06,14,45,00);
% titl = 'Mid Tide';
start = datenum(2015,03,06,15,20,00);     %high
stop = datenum(2015,03,06,15,30,00);
titl = 'High Tide';

%%%Data Structure%%%
dat = struct();

%%%Load ADCPS%%%
adcpdir = 'd:\Projects\Mekong_W2015\Data\Aquadopp\F2F2\';
toload = {'HR3_7March2015';'AD5116_9March2015';'AD5117_9March2015'};

for i = 1:3
    disp(['Loading ' toload{i}])
    load([adcpdir toload{i} '.mat'])
    ind = find(aqdp.datenum >= start & aqdp.datenum <= stop);
    w = aqdp.w(ind,:);wbar = nanmean(w,2); %depth average
    clear aqdp
    wbar = cmgbridge(wbar,10,100,1000);
    name = regexprep(toload{i},'(_[^_]*)','');
    
    %Spectra settings
    min_f = 1E-1;max_f = 1E2;
    n = length(wbar);
    fs = 8;
    nfft = floor(n/6);
    window = hanning(floor(n/10));
    noverlap = nfft/2;
    [pf,f] = pwelch(wbar,window,noverlap,nfft,fs);
    ff = min(find(f>=min_f));
    lf = max(find(f<=max_f));
    pff = pf(ff:lf);freq = f(ff:lf);
    dat.(name).f = freq;
    dat.(name).pf = pff;
end

%%%Load ADVs%%%
advdir = 'd:\Projects\Mekong_W2015\Data\Vector\F2F2\';
toload = {'V5109_050315';'VC01_050315';'VC101_050315'};

for i = 1:3
    disp(['Loading ' toload{i}])
    load([advdir toload{i} '.mat'])
    ind = find(ADV.datetime >= start & ADV.datetime <= stop);
    w = ADV.W(ind,:);wbar = nanmean(w,2); %depth average
    clear ADV
    wbar = cmgbridge(wbar,10,100,1000);
    name = regexprep(toload{i},'(_[^_]*)','');
    
    %Spectra settings
    min_f = 1E-1;max_f = 1E2;
    n = length(wbar);
    fs = 32;
    nfft = floor(n/6);
    window = hanning(floor(n/10));
    noverlap = nfft/2;
    [pf,f] = pwelch(wbar,window,noverlap,nfft,fs);
    ff = min(find(f>=min_f));
    lf = max(find(f<=max_f));
    pff = pf(ff:lf);freq = f(ff:lf);
    dat.(name).f = freq;
    dat.(name).pf = pff;
end

%%%Plot Routine%%%
fn = fieldnames(dat);
order = [1 2 4 5 3 6];
 
f1 = figure(1);
set(f1,'PaperOrientation','portrait',...
    'position',[400 50   500 1000]);
set(gcf,'color','w','PaperPositionMode','auto')
sp = zeros(1,6);
for i = 1:6
    sp(i) = subplot(6,1,i);
    f = dat.(fn{order(i)}).f;
    pf = dat.(fn{order(i)}).pf;
    loglog(f,pf,...
        'LineWidth',1.5,...
        'Color','k')
    text(4,1E-1,fn{order(i)},...
        'FontSize',13,'FontName','Cambria')
    ylabel('PSD','FontSize',13,'FontName','Cambria')
end
%Global Adjustments%
set(sp,'Xlim',[1E-1 1E1],...
    'Ylim',[1E-6 1E0],...
    'LineWidth',1.5,...
    'FontSize',13,...
    'FontName','Cambria',...
    'TickDir','out')
set(sp(1:5),'XTickLabel',[])
set(sp(1),'Position',[0.15 0.83 0.8 0.13])
set(sp(2),'Position',[0.15 0.68 0.8 0.13])
set(sp(3),'Position',[0.15 0.53 0.8 0.13])
set(sp(4),'Position',[0.15 0.38 0.8 0.13])
set(sp(5),'Position',[0.15 0.23 0.8 0.13])
set(sp(6),'Position',[0.15 0.08 0.8 0.13])
xlabel(sp(6),'Frequency (Hz)')
figdir = 'd:\Projects\Mekong_W2015\Figures\TransectSpectra\';
export_fig([figdir titl],'-jpeg','-nocrop')

