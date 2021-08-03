%Code to plot X-Wavelet Power Spectra, copied from wtc.m
%Instructions: Run wtc on data, save [Rsq,period,coi,wtcsig,t,Wxy], run
%this script.
function plotwtc(Rsq,period,coi,wtcsig,t,Wxy,dt)
%Default arguments for plot routine
Args=struct('ArrowDensity',[30 30],...
            'ArrowSize',1.5,...
            'ArrowHeadSize',1.5);
ad=mean(Args.ArrowDensity);
Args.ArrowSize=Args.ArrowSize*30*.03/ad;
%Args.ArrowHeadSize=Args.ArrowHeadSize*Args.ArrowSize*220;
Args.ArrowHeadSize=Args.ArrowHeadSize*120/ad;

Yticks = 2.^(fix(log2(min(period))):fix(log2(max(period))));

H=imagesc(t,log2(period),Rsq);
%[c,H]=safecontourf(t,log2(period),Rsq,[0:.05:1]);
%set(H,'linestyle','none')

set(gca,'clim',[0 1])

% HCB=colorbar;

set(gca,'YLim',log2([min(period),max(period)]), ...
    'YDir','reverse', 'layer','top', ...
    'YTick',log2(Yticks(:)), ...
    'YTickLabel',num2str(Yticks'), ...
    'layer','top')
% ylabel('Period')
hold on

%phase plot
aWxy=angle(Wxy);
aaa=aWxy;
aaa(Rsq<.5)=NaN; %remove phase indication where Rsq is low
%[xx,yy]=meshgrid(t(1:5:end),log2(period));

phs_dt=round(length(t)/Args.ArrowDensity(1)); tidx=max(floor(phs_dt/2),1):phs_dt:length(t);
phs_dp=round(length(period)/Args.ArrowDensity(2)); pidx=max(floor(phs_dp/2),1):phs_dp:length(period);
phaseplot(t(tidx),log2(period(pidx)),aaa(pidx,tidx),Args.ArrowSize,Args.ArrowHeadSize);

if ~all(isnan(wtcsig))
    [~,h] = contour(t,log2(period),wtcsig,[1 1],'k');
    set(h,'linewidth',2)
end
%suptitle([sTitle ' coherence']);
tt=[t([1 1])-dt*.5;t;t([end end])+dt*.5];
hcoi=fill(tt,log2([period([end 1]) coi period([1 end])]),'w');
% set(hcoi,'alphadatamapping','direct','facealpha',.5)
hold off
end