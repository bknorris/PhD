function plotaqdpvels(aqdp,ornt,hcaxis,vcaxis)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Plot aquadopp velocities, correlations and backscatter
%
%  Usage:
%  plotaqdpvels(aqdp,ornt,hcaxis,vcaxis)
%
%  Where ornt = instrument orientation (UP or DOWN)
%        hcaxis = horizontal color axis specification for velocity color
%        scaling (ex. 0-1)
%        vcaxis = verticalcolor axis specification for velocity color
%        scaling
%
%  Script developed by Benjamin K Norris, University of Waikato, NZ 2014
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
yearday = aqdp.yearday;
rangebins = aqdp.rangebins;
theFillValue = 0;
date = datestr(aqdp.datenum(1),'ddmmyyyy');
hrs = (yearday-(floor(yearday(1))))*24;
corav = (aqdp.cor1+aqdp.cor2+aqdp.cor3)/3;
backav = (aqdp.beam1+aqdp.beam2+aqdp.beam3)/3;
ind = isnan(aqdp.u);aqdp.u(ind) = theFillValue;u = aqdp.u;
ind = isnan(aqdp.v);aqdp.v(ind) = theFillValue;v = aqdp.v;
ind = isnan(aqdp.w);aqdp.w(ind) = theFillValue;w = aqdp.w;
ind = isnan(corav);corav(ind) = theFillValue;
ind = isnan(backav);backav(ind) = theFillValue;


max1 = hcaxis;
min1 = -hcaxis;
step1 = (-min1+max1)/4;
max2 = vcaxis;
min2 = -vcaxis;
step2 = (-min2+max2)/4;

f1 = figure(1);
set(f1,'PaperOrientation','portrait',...
    'position',[200 200   1250   550]);
a = tight_subplot(5,1,0.05,[0.1 0.05],[0.1 0.05]);

axes(a(1))
if strcmp(ornt,'DOWN');
    imagesc(hrs,flipud(rangebins),flipud(u'))
    hold on
    p(1) = plot(hrs,aqdp.pressure,'k','LineWidth',1);
else
    imagesc(hrs,rangebins,u')
    hold on
    p(1) = plot(hrs,aqdp.pressure,'k','LineWidth',1);
end
title('\bf\itU - East Velocity')
caxis([min1 max1])
c = colorbar;
cbh = get(c,'Title');
titleString = '[m/s]';
set(cbh ,'String',titleString);
set(c,'YTick',[min1:step1:max1]);
g = gca;
set(g,'XTickLabel',[])
ylabel('\bf\itrange (m)')

axes(a(2))
if strcmp(ornt,'DOWN');
    imagesc(hrs,flipud(rangebins),flipud(v'))
    hold on
    p(2) = plot(hrs,aqdp.pressure,'k','LineWidth',1);
else
    imagesc(hrs,rangebins,v')
    hold on
    p(2) = plot(hrs,aqdp.pressure,'k','LineWidth',1);
end
title('\bf\itV - North Velocity')
caxis([min1 max1])
c = colorbar;
cbh = get(c,'Title');
titleString = '[m/s]';
set(cbh ,'String',titleString);
set(c,'YTick',[min1:step1:max1]);
g = gca;
set(g,'XTickLabel',[])
ylabel('\bf\itrange (m)')

axes(a(3))
if strcmp(ornt,'DOWN');
    imagesc(hrs,flipud(rangebins),flipud(w'))
    hold on
    p(3) = plot(hrs,aqdp.pressure,'k','LineWidth',1);
else
    imagesc(hrs,rangebins,w')
    hold on
    p(3) = plot(hrs,aqdp.pressure,'k','LineWidth',1);
end
title('\bf\itW - Vertical Velocity')
caxis([min2 max2])
c = colorbar;
cbh = get(c,'Title');
titleString = '[m/s]';
set(cbh ,'String',titleString);
set(c,'YTick',[min2:step2:max2]);
g = gca;
set(g,'XTickLabel',[])
ylabel('\bf\itrange (m)')

axes(a(4))
if strcmp(ornt,'DOWN');
    imagesc(hrs,flipud(rangebins),flipud(corav'))
    hold on
    p(4) = plot(hrs,aqdp.pressure,'k','LineWidth',1);
else
    imagesc(hrs,rangebins,corav')
    hold on
    p(4) = plot(hrs,aqdp.pressure,'k','LineWidth',1);
end
title('\bf\itAverage Correlations')
caxis([0 100])
c = colorbar;
cbh = get(c,'Title');
titleString = '[%]';
set(cbh ,'String',titleString);
set(c,'YTick',[0:25:100]);
g = gca;
set(g,'XTickLabel',[])
ylabel('\bf\it%')

axes(a(5))
if strcmp(ornt,'DOWN');
    imagesc(hrs,flipud(rangebins),flipud(backav'))
    hold on
    p(5) = plot(hrs,aqdp.pressure,'k','LineWidth',1);
else
    imagesc(hrs,rangebins,backav')
    hold on
    p(5) = plot(hrs,aqdp.pressure,'k','LineWidth',1);
end
title('\bf\itAverage Backscatter')
c = colorbar;
cbh = get(c,'Title');
titleString = '[Counts]';
set(cbh ,'String',titleString);
ylabel('\bf\itCounts')
xlabel(['\bf\itHours Elapsed from ' date])
set(a,...
    'YTick',[0:0.1:0.5],...
    'YLim', [0 0.5],...
    'XGrid','on',...
    'GridLineStyle','-',...
    'YDir','normal')
set(gcf,'color','w','PaperPositionMode','auto');

filen = regexprep(aqdp.metadata.name,'_f','');
fdir = 'c:\Users\bkn5\Projects\';
fpath = [pwd, '\',mfilename];whichdir = strfind(fpath,'Mekong');
folder = fpath(whichdir:whichdir+11);
fpath = [fdir folder '\Figures\'];
if ~exist([fpath filen],'dir')
    [a,~,~] = mkdir(fpath,filen);
    if a == 1
        disp(['Figure Directory ' fpath filen '\' ' created'])
    end
end
filep = [fpath filen '\'];
fname = [filen '_' datestr(aqdp.datenum(1),'ddmmyy')];
if exist([filep fname '.png'],'file')
    pause(4)
    disp(['Figure ' fname '.png exists at this location'])
    prompt = 'Overwrite figure? [y/n] ';
    result = input(prompt,'s');
    if strcmp(result,'y') || strcmp(result,'yes');
        export_fig([filep fname],'-png','-m1','-r900','-opengl')
        disp(['Figure ' fname '.png saved'])
    end
    if strcmp(result,'n') || strcmp(result,'no');
        disp('Figure not saved')
    end
else
    export_fig([filep fname],'-png','-m1','-r900','-opengl')
    disp(['Figure ' fname '.png saved'])
end
end


