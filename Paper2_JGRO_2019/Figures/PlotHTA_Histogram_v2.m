%Correlates peaks in turbulence with peaks in rotated (Vectrino) velocities 
%and plots on a histogram for onshore/offshore flow. I've left the Vector
%load script in to plot wave statistics.

clear
%load vectrino/turbulence files from HTA1
load('d:\Projects\Mekong_W2015\DataAnalysis\Paper2\Turbulence\7March2015_TKEbdadj.mat')
load('d:\Projects\Mekong_W2015\DataAnalysis\Paper2\AveragedVelocities\7March2015_bdadj.mat')
fn = fieldnames(Stat);
data = struct();
%%%Find peaks in Turbulence T-S%%%
pkhght = 1E-4; %minimum peak in TKE
pkdist = 2; %peaks must be > 240 sec apart
for i = 1:3
    E = (Stat.(fn{i}).z1.E+Stat.(fn{i}).z2.E)./2;
    eps = nanmean(E(1:5,:));
    [pks,lc] = findpeaks(eps,...
        'minpeakheight',pkhght,...
        'minpeakdistance',pkdist);
    pkvels = Avgs.(fn{i}).x(lc);
    data.(fn{i}).pkvels =  pkvels;
    data.(fn{i}).peaks = pks;
    data.(fn{i}).loc = lc;
    data.(fn{i}).eps = eps;
end

%%%Plot Figures%%%
f1 = figure(1);
set(f1,'PaperOrientation','portrait',...
    'position',[400 200   500  400]);
set(gcf,'color','w','paperpositionmode','auto')
c = flipud([0.1 0.1 0.1;0.5 0.5 0.5;0.7 0.7 0.7]);
ct = [1 4 7];
for i = 1:3
    cts = ct(i);
    lc = data.(fn{i}).loc;
    pks = data.(fn{i}).peaks;
    pkvels = data.(fn{i}).pkvels;
    s = sign(pkvels);
    neg = sum(s(:)==-1);
    pos = sum(s(:)==1);
    center = [cts-0.5 cts+0.5];
    br = bar(center,[neg pos]); hold on
    set(br,'facecolor',c(i,:),...
        'linewidth',1.5)
    text(cts-0.7,240,'Off     On')
    %UPDATE: 06/04/2017 calculate mean(eps) of peaks
    disp(fn{i})
    eps = nanmean(data.(fn{i}).eps(s<0));
    disp(['mean offshore eps: ' sprintf('%.2d',eps) ' W/kg'])
    eps = nanmean(data.(fn{i}).eps(s>0));
    disp(['mean onshore eps: ' sprintf('%.2d',eps) ' W/kg'])
end
break
set(gca,'ylim',[0 400],'ytick',0:100:400,...
    'xticklabel',{'','x = -10cm','','',...
    'x = 10cm','','','x = 20cm'})
ylabel('No. of Occurences')
%positioning
set(sp(1),'position',[0.08 0.56 0.42 0.3])
set(sp(2),'position',[0.08 0.15 0.42 0.3])
set(sp(3),'position',[0.575 0.1 0.4 0.83])
prettyfigures('text',13,'labels',14,'box',1)
set(leg,'position',[0.11 0.34 0.05 0.05],...
    'box','off')
savefigdir = 'd:\Projects\Mekong_W2015\Figures\Paper2\WaveVelsTurbulence\';
% export_fig([savefigdir 'Vels_eps_hist_v2'],'-pdf','-nocrop')