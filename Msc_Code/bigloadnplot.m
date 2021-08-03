%Big loadnplot for all the SAND runs

clear
filelist = dir('*PAR.mat');
files = {filelist.name};
dirname = 'D:\Projects\PARSensor\Clay';
%Do the load thing
names = {'blk1000';'blk100';'blk10';'red1000';'red100';'red10';'tan1000';'tan100';'tan10';...
    'wht1000';'wht100';'wht10'}';
data = {};
for i = 1:length(files)
    fname = fullfile(dirname,files{i});
    data{i} = load(fname);
end
clay = cell2struct(data,names,2);
depth = 0:5:160; %step interval (in cm)

%one quick issue with red10:
clay.red10.par.average = clay.red10.par.average*-1;
%plotting classifications:
tens = [3 6 9 12];
onehunned = [2 5 8 11];
mil = [1 4 7 10];

cmap = [0 0 0; 0 0 0; 0 0 0; 1 0 0; 1 0 0; 1 0 0; 0.6 0.4 0.2; 0.6 0.4 0.2; 0.6 0.4 0.2;...
    1 1 1; 1 1 1; 1 1 1];

m1 = figure('PaperOrientation','portrait',...
    'position',[400   100   800   800]);

for i=mil
    t(i) = plot(clay.(names{i}).par.average,depth,...
        'color',cmap(i,:),...
        'LineWidth',2,...
        'MarkerSize',8);
    hold on
end
t(t==0) = [];
leg = legend(t,names(mil),'location','southeast');
set(leg,'box','off')

haxes1 = gca;
set(haxes1,...
    'YDir','Reverse',...
    'Color',[0.8 0.8 0.8],...
    'YGrid','on',...
    'TickDir','in',...
    'XAxisLocation','bottom',...
    'XColor',[.2 .2 .2], ...
    'YColor',[.2 .2 .2], ...
    'LineWidth',1,...
    'TickLength',[.02 .02],...
    'YMinorTick','on',...
    'YTick',[0:10:160],...
    'YLim', [0 160])
title('\bfPAR Attenuation by Sediment Color, 1000mg/L Assumed Concentration')
ylabel(haxes1,'\itDepth (cm)')
xlabel(haxes1,'\itW/m^2')
box on
set(gcf, 'PaperPositionMode', 'auto');
hgexport(gcf, '1000mg.jpg', hgexport('factorystyle'), 'Format', 'jpeg');

    

    