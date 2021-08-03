%Generate Bed Trace for VPs for spectral analysis
%
% Paper 3: Sediment Motion in Mangroves
%
% This is Version 1.0 of this script
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear
close all

%let's start with F2F2_2 as an example
dir1 = 'e:\Mekong_W2015\DataAnalysis\Paper3\VPs\06-03-15\';
files = {'6March2015_Vels.mat';'6March2015_Sen.mat'};
hab = 0.063;
expname = 'F2F2_2';
fn1 = whos('-file',[dir1 files{1}]);fn1 = {fn1.name};
V = matfile([dir1 files{1}]);
S = matfile([dir1 files{2}]);
% load('d:\Projects\Mekong_W2015\Data\Aquadopp\F2F2\AD5116_9March2015.mat')

%load the data
for ii = 1
    vp = V.(fn1{ii});
    se = S.(fn1{ii});
    
    bdt = vp.bdtime;
    bdist = vp.bdist;
    amp = (se.Amp1+se.Amp2+se.Amp3+se.Amp4)./4;
    ampt = se.time;
    vph = hab-vp.rb;
    figure(1) %plot amp
    subplot(121)
    imagesc(ampt,vph,amp'),set(gca,'ydir','normal')
    hold on
    plot(bdt,hab-bdist)
    datetick('x','HH:MM','keepticks','keeplimits')
    xlabel(['Time on ' datestr(bdt(1),'dd-mm-yy')])
    ylabel('Height Above Bottom (m)')
    set(gca,'ylim',[min(hab-bdist) hab])
    
    sfiledir = 'g:\Mekong_W2015\DataAnalysis\Paper3\BottomTrack\06-03-15\';
    bds = my_running_max(bdist,50);
    bd = my_running_median(bds,1024);
    subplot(122)
    imagesc(ampt,vph,amp'),set(gca,'ydir','normal')
    hold on
    plot(bdt,hab-bd)
    datetick('x','HH:MM','keepticks','keeplimits')
    xlabel(['Time on ' datestr(bdt(1),'dd-mm-yy')])
    set(gca,'ylim',[min(hab-bdist) hab])
    sfdir = 'e:\MemoryStick\GradSchool\DataAnalysis\Paper3\Figures\';
    prettyfigures('text',12,'labels',13,'box',1,'gcolor','k')
    export_fig([sfdir 'Backscatter_BLE_060315'],'-pdf','-nocrop')

    btrace.(fn1{ii}).time = bdt;
    btrace.(fn1{ii}).bdist = bd;

end
% save([sfiledir expname '_bdtrace_v2'],'-struct','btrace')