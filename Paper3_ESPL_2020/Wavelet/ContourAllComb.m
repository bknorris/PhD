cols = {'eventl';'sigh';'depth';'orbwv';'umag';'usqd';'taub';'eps'};
cids = 1:length(cols);
cmbid = combvec(cids,cids);
for ii = 1:length(cmbid)
    if cmbid(1,ii) == cmbid(2,ii)
        continue
    else
        xvar = cols{cmbid(1,ii)};
        yvar = cols{cmbid(2,ii)};
    end
    %%%
    sp = zeros(3,1);
    f = figure(ii);
    set(f,'PaperOrientation','portrait',...
        'position',[400 100   1000   400],...
        'renderer','painters');
    factor = 20; %joss's npred
    range = 1; %joss's perc_range
    cmap = brewermap(factor,'RdBu');cmap = flipud(cmap);
    for i = 1:3
        sp(i) = subplot(1,3,i);
        %create a patch for background color
        %     px = [0 1.5 1.5 0];
        %     py = [0 0 1.5 1.5];
        %     patch(px,py,'k'),hold on
        xb = [dat.(fn{i}).wave.(xvar); dat.(fn{i}).ig.(xvar)];
        yb = [dat.(fn{i}).wave.(yvar); dat.(fn{i}).ig.(yvar)];
        bdmed = [dat.(fn{i}).wave.bdmed; dat.(fn{i}).ig.bdmed];
        %prep data for contour
        mx = linspace(min(xb)*range,max(xb)*range,factor);
        my = linspace(min(yb)*range,max(yb)*range,factor);
        [vx,vy] = meshgrid(mx,my);
        v = griddata(xb,yb,bdmed,vx,vy,'cubic');
        p = contourf(vx,vy,v);
        colormap(cmap);
        % shading(gca,'interp')
        caxis([-0.03 0.03])
    end
    title([xvar ' vs. ' yvar])
end