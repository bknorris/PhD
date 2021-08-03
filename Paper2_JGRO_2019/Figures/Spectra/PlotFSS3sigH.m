%Plot wave heights from FSS3, mudflat and frige.
%This is to show that there is little wave height descrepancy between the 
%unvegetated mudflat and the vegetated fringe environment (i.e. importance 
%of vegetation to TKE measurements and WPF reduction). 

clear
ddir = '/Volumes/MEKONG1/Mekong_W2015/DataAnalysis/WavePowerFlux/TestBed/';
DATA.mudflat1 = load([ddir 'WpfV5108_FSS31.mat']);
DATA.fringe1 = load([ddir 'WpfAD5116_FSS31.mat']);
DATA.forest1 = load([ddir 'WpfVC101_FSS31.mat']);

DATA.mudflat2 = load([ddir 'WpfADHR3_FSS32.mat']);
DATA.fringe2 = load([ddir 'WpfAD5116_FSS32.mat']);
DATA.forest2 = load([ddir 'WpfVC101_FSS32.mat']);

ddir = '/Volumes/MEKONG1/Mekong_W2015/DataAnalysis/WavePowerFlux/FSS3.3/';
DATA.mudflat3 = load([ddir 'WpfADHR3_FSS33.mat']);
DATA.fringe3 = load([ddir 'WpfAD5116_FSS33.mat']);
DATA.forest3 = load([ddir 'WpfVC101_FSS33.mat']);

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
set(gcf,'color','w','PaperPositionMode','auto')
sp(1) = subplot(311); hold on
p(1) = plot(Waves.mudflat1.time,Waves.mudflat1.Hsig,'-kd');
p(2) = plot(Waves.fringe1.time,Waves.fringe1.Hsig,'--ro');
p(3) = plot(Waves.forest1.time,Waves.forest1.Hsig,'-.b^');
set(gca,'XTickLabel',[])
ylabel('\bf\itH sig (m)','FontSize',16,'FontName','Helvetica')
text(0.5,0.9,'\bf\itDay 1','units','normalized','FontSize',16)

sp(2) = subplot(312); hold on
p(4) = plot(Waves.mudflat2.time,Waves.mudflat2.Hsig,'-kd');
p(5) = plot(Waves.fringe2.time,Waves.fringe2.Hsig,'--ro');
p(6) = plot(Waves.forest2.time,Waves.forest2.Hsig,'-.b^');
set(gca,'XTickLabel',[])
ylabel('\bf\itH sig (m)','FontSize',16,'FontName','Helvetica')
text(0.5,0.9,'\bf\itDay 2','units','normalized','FontSize',16)

sp(3) = subplot(313); hold on
p(7) = plot(Waves.mudflat3.time,Waves.mudflat3.Hsig,'-kd');
p(8) = plot(Waves.fringe3.time,Waves.fringe3.Hsig,'--ro');
p(9) = plot(Waves.forest3.time,Waves.forest3.Hsig,'-.b^');
ylabel('\bf\itH sig (m)','FontSize',16,'FontName','Helvetica')
xlabel('\bf\itMinutes After Low Tide','FontSize',16,'FontName','Helvetica')
text(0.5,0.9,'\bf\itDay 3','units','normalized','FontSize',16)
set(p,'LineWidth',1.5,'MarkerSize',8,'MarkerFaceColor',[1 1 1]);
%global plot adjustments
set(sp,'LineWidth',1.5,'FontSize',16,'FontName','Helvetica',...
    'XTick',0:ts*2:length(Waves.mudflat2.Hsig)*ts,...
    'Xlim',[0 length(Waves.mudflat2.Hsig)*ts],...
    'Ylim',[0 0.4],'box','on')
set(sp(1),'position',[0.125 0.63 0.79 0.25])
set(sp(2),'position',[0.125 0.355 0.79 0.25])
set(sp(3),'position',[0.125 0.08 0.79 0.25])

leg = legend([p(1) p(2) p(3)],{'\bf\itMudflat','\bf\itFringe','\bf\itForest'},...
    'box','on');
set(leg,'position',[0.79 0.45 0.15 0.1],'LineWidth',1.5)

fpath = '/Volumes/MEKONG1/Mekong_W2015/Figures/Spectra&Waves/';
fname = 'HsigFSS3';
export_fig([fpath fname],'-pdf','-nocrop','-m1','-r900')
disp(['Figure ' fname '.pdf saved'])

