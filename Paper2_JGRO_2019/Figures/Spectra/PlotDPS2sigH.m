%Plot wave heights from DPS2, mudflat, fringe and forest. This is to show 
%that there is little wave height descrepancy between the unvegetated 
%mudflat and the vegetated fringe environment (i.e. importance of vegetation to TKE
%measurements and WPF reduction). 

clear
ddir = '/Volumes/MEKONG1/Mekong_W2015/DataAnalysis/WavePowerFlux/DPS2Turbulence/';
DATA.mudflat1 = load([ddir 'WpfV5109_DPS2.mat']);
DATA.fringe1 = load([ddir 'WpfAD5116_DPS2.mat']);
DATA.forest1 = load([ddir 'WpfAD5117_DPS2.mat']);

%timestep is 10 minutes
ts = 10;
sname = fieldnames(DATA);
Waves = struct();
for i = 1:length(sname);
    Hs = max(DATA.(sname{i}).WPF.Hs,[],2);
    idx = ts:ts:length(Hs);
    idx = [1 idx];
    for j = 1:length(idx)-1
        Hsig(:,j) = mean(Hs(idx(j):idx(j+1)));
    end
    Waves.(sname{i}).Hsig = Hsig;
    Waves.(sname{i}).time = idx(2:end);
    clear idx Hsig Hs
end

f1 = figure;    
set(f1,'PaperOrientation','portrait',...
    'position',[500 60   800 1050]);
set(gcf,'color','w','PaperPositionMode','auto'), hold on
p(1) = plot(Waves.mudflat1.time,Waves.mudflat1.Hsig,'-kd');
p(2) = plot(Waves.fringe1.time,Waves.fringe1.Hsig,'--ro');
p(3) = plot(Waves.forest1.time,Waves.forest1.Hsig,'-.b^');
ylabel('\bf\itH sig (m)','FontSize',16,'FontName','Helvetica')
xlabel('\bf\itMinutes After Low Tide','FontSize',16,'FontName','Helvetica')
set(p,'LineWidth',1.5,'MarkerSize',8,'MarkerFaceColor',[1 1 1]);
%global plot adjustments
set(gca,'LineWidth',1.5,'FontSize',16,'FontName','Helvetica',...
    'XTick',0:ts*2:length(Waves.mudflat1.Hsig)*ts,...
    'Xlim',[0 length(Waves.mudflat1.Hsig)*ts],...
    'Ylim',[0.04 0.1],'box','on')

leg = legend([p(1) p(2) p(3)],{'\bf\itMudflat','\bf\itFringe','\bf\itForest'},...
    'box','on');
set(leg,'position',[0.72 0.15 0.15 0.1],'LineWidth',1.5)

fpath = '/Volumes/MEKONG1/Mekong_W2015/Figures/Spectra&Waves/';
fname = 'HsigDPS2';
export_fig([fpath fname],'-pdf','-nocrop','-m1','-r900')
disp(['Figure ' fname '.pdf saved'])

