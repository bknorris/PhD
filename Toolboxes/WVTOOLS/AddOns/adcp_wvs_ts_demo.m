% adcp_wvs_ts_demo.m  An mfile to plot the timeseries of wave parameters from
%                     ADCP data.
%
% Copyright 2004 
% USGS Woods Hole Field Center
% Written by Charlene Sullivan
% csullivan@usgs.gov

% C. Sullivan 10/21/04, Raw ADCP data is first processed through WavesMon
%                       software.  The mfile adcpWvs2nc.m transforms the
%                       data output by WavesMon into a netCDF file.  It is
%                       assumed that the netCDF files are stored in a 
%                       directory named according to mooring number and 
%                       instrument type.
%             02/07/05, Datafiles now stored in 'dataDirectory'.  
%                       Interpolate data onto a standard time axis.  Name
%                       pdf's waves_site#.pdf. 
%             03/03/05, Reformatting of labels and ticks
%             07/29/05, This mfile is only a DEMO version and is coded to
%                       work with a specific data set.
% Dependencies:
%       the netcdf toolbox

dataDirectory='C:\waves_module\722\7221wh\wavesmon';
plotDirectory='C:\waves_module\722\7221wh\wavesmon';

% define data filenames 
dataFiles='7221adwp-cal.nc';
       
% define parameters for plots
col=get(0,'defaultaxescolororder');                 %plot color order
ylabels={ {'Significant','Wave Height (m)'},...                    %plot ylabels
          {'Peak Wave','Period (s)'},...
          {'Wave','Direction (deg)'},...
          {'Water','Depth (m)'},...
          {'Frequency (Hz)'} };
startTime=julian(2003,10,28,00);                            %plot start stime
stopTime=julian(2003,11,17,00);                             %plot stop time  
timediff=stopTime-startTime;
xlims=[startTime stopTime];
xt=[startTime:5:stopTime];

% initialize variables
jd=cell(length(dataFiles),1);
hs=cell(length(dataFiles),1);
tp=cell(length(dataFiles),1);
wvdir=cell(length(dataFiles),1);
depth=cell(length(dataFiles),1);
vspec=cell(length(dataFiles),1);

% load data
nc=netcdf(dataFiles, 'nowrite');
jd=nc{'time'}(:) + (nc{'time2'}(:)/3600/1000/24);
hs=nc{'wh_4061'}(:);                     %significant wave height, meters
hs(hs>1e34)=nan;
tp=nc{'wp_peak'}(:);                     %peak period, seconds
tp(tp>1e34)=nan;
wvdir=nc{'wvdir'}(:);                    %peak wave direction, degrees
wvdir(wvdir>1e34)=nan;
depth=nc{'hght_18'}(:);                  %water depth, meters
depth(depth>1e34)=nan;
vspec=nc{'vspec'}(:)./1000;              %non-directional velocity-derived spectra, M/sqrt(Hertz)
vspec(vspec==-32768/1000)=nan;
freq=nc{'frequency'}(:);                         %frequency, Hertz
site=nc.DESCRIPT(6);                        %site number
theFillValue=nc{'u_1205'}.FillValue_(:);    %netcdf fill value
delta_t = gmean(diff(jd));               %time step, julian days
nc=close(nc);

% plot data
figure('position',[100 100 1000 800])
han(1)=subplot(5,1,1);
plot(jd,hs)
ylim([0 2])
title(['Site ',num2str(site),' Wave Observations'],'fontsize',14)
han(2)=subplot(5,1,2);
plot(jd,tp)
ylim([0 20])
han(3)=subplot(5,1,3);
plot(jd,wvdir)
ylim([0 360])
han(4)=subplot(5,1,4);
plot(jd,depth)
ylim([7.5 10.5])
han(5)=subplot(5,1,5);
jd=repmat(jd',length(freq),1);
freq=repmat(freq,1,length(jd(1,:)));
pcolor(jd,freq,vspec')
set(gca,'Layer','top')
ylim([0 .5])
shading flat
caxis([-0.5 1.5])
axpos=get(han(5),'Position');
cbar=colorbar;
cbarpos=get(cbar,'Position');
set(cbar,'Position',[0.91 axpos(2) cbarpos(4)/4 cbarpos(4)])
set(get(cbar,'Ylabel'),'String','M/sqrt(Hz)')

% reformat xlims, and ylabels
for h=1:length(han)
    axes(han(h))
    set(gca,'xlim',xlims,'xtick',xt,'xticklabel',[])
    ylabhan=get(gca,'Ylabel');
    set(ylabhan,'String',ylabels{h})
    pos=get(ylabhan,'position');
    set(get(gca,'Ylabel'),'position',[julian(2003,10,27) pos(2) pos(3)])
end

% re-do xlabels
axes(han(5))
oldhan=findobj(gca,'Tag','xtick labels');
delete(oldhan);
for tic=1:length(xt)
    xtg=gregorian(xt(tic));
    xtd=datenum(xtg);
    xtl={[datestr(xtd,3),' ',datestr(xtd,7)];datestr(xtd,10)};
    xtlpos=ylims(1,1)-0.1;
    text(xt(tic),xtlpos,xtl,'horizontalalignment','center',...
        'verticalalignment','top')
end

% the ouput pdf filename
orient landscape
set(gcf,'PaperPositionMode','Auto')
set(gcf,'Renderer','Painters')
print(fullfile(plotDirectory,['waves_site',site]),'-dpdf')


