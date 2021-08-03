function rss=reynolds_stress( au,av,aw,bu,bv,bw,fs,bn,jd,iverbose )
% REYNOLDS_STRESS - Calculate <u'w'> and <v'w'> using Trowbridge method
% rss=reynolds_stress(au,av,aw,bu,bv,bw,fs,bn,jd,iverbose)
%
% Returns:
%       rss.rsu
%       rss.rsv
%       rss.mag
%       rss.dir
%       rss.e2x
%       rss.e2y
%       rss.e2
%       rss.dof
% Rotation does not matter for the bulk estimate, but rotation is used in
% the error estimates
global figc
if(exist('iverbose')~=1),iverbose=0,end;
DirA = 0;     %not sure what this was used for in calling program
ialign = 0;
zref = 0.35; %reference height for measurements (used in calc. dof)
% simplest estimate
xd = au(:)-bu(:);
yd = av(:)-bv(:);
zd = aw(:)-bw(:);
cxw = cov([xd zd],1);
rsu = -0.5*cxw(1,2);
cyw = cov([yd zd],1);
rsv = -0.5*cyw(1,2);
% return values
rss.rsu = rsu;
rss.rsv = rsv;
[rs_mag,rs_dir]=pcoord(rsu,rsv);
rss.mag = rs_mag;
rss.dir = rs_dir;

[sam dam]=pcoord(mean(au),mean(av));
[sbm dbm]=pcoord(mean(bu),mean(bv));
[sabm dabm]=pcoord(mean([au(:);bu(:)]),mean([av(:);bv(:)]));
% rotate so u = mean directio of A
[sa da]=pcoord(au(:),av(:));
[sb db]=pcoord(bu(:),bv(:));
% rotate to +u = mean flow direction    db_r = db-DirA+90;
da_r = da-dabm+90;
db_r = db-dabm+90;
[u_ar,v_ar]=xycoord(sa,da_r);
[u_br,v_br]=xycoord(sb,db_r);
[sarm darm]=pcoord(mean(u_ar),mean(v_ar));
[sbrm dbrm]=pcoord(mean(u_br),mean(v_br));
xd = u_ar(:)-u_br(:);
yd = v_ar(:)-v_br(:);
zd = aw(:)-bw(:);
% estimates in alongstream coordinates
cxw = cov([xd zd],1);
rsur = -0.5*cxw(1,2);
cyw = cov([yd zd],1);
rsvr = -0.5*cyw(1,2);

dof = sabm*(length(xd))./(fs*zref);
e1x = 0.5*( var(xd)*var(zd) )/sqrt(dof);
e1y = 0.5*( var(yd)*var(zd) )/sqrt(dof);
e1 = sqrt( e1x.^2 + e1y.^2);
% This is what I get from Bendat & Piersol, eqn. 8.103
e2x = rsur*0.5*sqrt( (1+(cxw(1,1)*cxw(2,2))./cxw(1,2).^2)./dof );
e2y = rsvr*0.5*sqrt( (1+(cyw(1,1)*cyw(2,2))./cyw(1,2).^2)./dof );
% This is more like what John T and Soulsby, eqn. 18 indicate
% e2x = rsur*0.5*sqrt( ((cxw(1,1)*cxw(2,2))./cxw(1,2).^2)./dof );
% e2y = rsvr*0.5*sqrt( ((cyw(1,1)*cyw(2,2))./cyw(1,2).^2)./dof );
e2 = sqrt( e2x.^2 + e2y.^2);
rss.e2x = e2x;
rss.e2y = e2y;
rss.e2 = e2;
rss.dof = dof;

if(iverbose)
   fprintf(1,'   sam, dam: %g %g\n',sam,dam)
   fprintf(1,'   sbm, dbm: %g %g\n',sbm,dbm)
   fprintf(1,' sabm, dabm: %g %g\n',sabm,dabm)
   fprintf(1,' sarm, darm: %g %g\n',sarm,darm)
   fprintf(1,'sbrbm, dbrm: %g %g\n',sbrm,dbrm)
   fprintf(1,' rsu,   rsv: %g %g\n',rsu,rsv);
   fprintf(1,'rsur, rsvr,: %g %g\n',rsur,rsvr)
   %   [rs_mag,rs_dir]=pcoord(rsu,rsv); % calc'd above
   fprintf(1,'  rs_mag, rs_dir: %g %g\n',rs_mag,rs_dir);
   [rsr_mag,rsr_dir]=pcoord(rsur,rsvr);
   fprintf(1,'rsr_mag, rsr_dir: %g %g\n',rsr_mag,rsr_dir);
   fprintf(1,'Trowbridge e       x, y, both: %g %g %g\n',e1x,e1y,e1);
   fprintf(1,'Bendat & Piersol e x, y, both: %g %g %g\n',e2x,e2y,e2);
   fprintf(1,'dof %5.1f\n',dof);
end
pct = -1;
if(iverbose) % do extra stuff
   n = length(xd);
   n4 = floor(n/4);
   [xxz,lags]=xcorr(xd,zd,n4,'coeff');
   [xyz,lags]=xcorr(yd,zd,n4,'coeff');
   [ax,lags]=xcorr(xd,n4,'coeff');
   [ay,lags]=xcorr(yd,n4,'coeff');
   [az,lags]=xcorr(zd,n4,'coeff');
   lag0 = find(lags==0);
   
   % this is not the time scale, but is consistent with dof
   timescale = zref*fs./sabm;
   fprintf(1,'time scale (lags) = %g\n',timescale)
   
   
   figure(2); clf
   subplot(211);
   h1=plot(lags(lag0:end)/fs,az(lag0:end),'-b','linewidth',1.5);
   hold on
   h2=plot(lags(1:lag0)/fs,ax(1:lag0),'-r','linewidth',1.5);
   h3=plot(lags/fs,xxz,'-k','linewidth',1.5);
   %h2=plot(lags/fs,xyz,'-c','linewidth',1.5);
   ts = sprintf('Burst %d, %s',bn,datestr(jd),0); % jd is actually dn
   %ts = sprintf('Burst %d, %s',bn,datestr(datenum(gregorian(jd)),0));
   text(100,.45,ts)
   ts = sprintf('Speed, Dir: %4.2f cm/s, %5.1f\\circT',100*sabm,dabm);
   text(100,.35,ts)
   ts = sprintf('{\\itu}_{*EC}: %5.2f+/-%5.2f cm/s',100*sqrt(rs_mag),100*sqrt(e2));
   text(100,.25,ts)
   Tx = cumsum(ax(lag0:end,1))./ax(lag0);
   Ty = cumsum(ay(lag0:end,1))./ay(lag0);
   Tw = cumsum(az(lag0:end,1))./az(lag0);
   legend([h1;h2;h3],'C_{zz}','C_{xx}','C_{xz}',...
      'location','northwest')
   axis([-300 300 -.2 .5])
   xlabel('Lag (s)')
   ylabel('Correlation')
   
   % COSPECPLOT -Plot cumulative cospectra
   nseg = 16;
   nfft = length(xd)/nseg;
   [Pxz,ff]=csd(xd,zd,nfft,fs,hanning(nfft),'mean');
   [Pyz,ff]=csd(yd,zd,nfft,fs,hanning(nfft),'mean');
   % TEST TEST TEST of replacement for csd:
   [Pxzp,ffp]=cpsd(xd,zd,hanning(nfft),[],nfft,fs);
   [Pyzp,ffp]=cpsd(yd,zd,hanning(nfft),[],nfft,fs);
   
   subplot(234)
   semilogx(ff,abs(Pxz),'-r','linewidth',2)
   hold on
   semilogx(ff,abs(Pyz),'-b','linewidth',2)
   semilogx(ffp,abs(Pxzp),'--m','linewidth',2)
   semilogx(ffp,abs(Pyzp),'--c','linewidth',2)
   plot([.02;10],[0 0],'--k')
   %axis([.02 10 -.05 .05])
   set(gca,'xlim',[.02 10])
   set(gca,'xtick',[.02 .1 1 10]);
   axis('square')
   % unsmoothed estimates for cumulative plots
   nseg = 1;
   nfft = length(xd)/nseg;
   [Pxz,ff]=csd(xd,zd,nfft,fs,boxcar(nfft),'mean');
   [Pyz,ff]=csd(yd,zd,nfft,fs,boxcar(nfft),'mean');
   [Pxzp,ffp]=cpsd(xd,zd,boxcar(nfft),0,nfft,fs);
   [Pyzp,ffp]=cpsd(yd,zd,boxcar(nfft),0,nfft,fs);
   % why is this 2 * delta f? CRS 1/5/2012.
   % dff = ff(3)-ff(1);
   % Instead, use:
   dff = ff(3)-ff(2);
   dffp = ffp(3)-ffp(2);
   ylabel('Co-spectra (m^2/s^2/Hz)')
   
   subplot(235)
   semilogx(ff,cumsum(abs(Pxz)*dff),'-r','linewidth',2)
   hold on
   semilogx(ff,cumsum(abs(Pyz)*dff),'-b','linewidth',2)
   semilogx(ffp,cumsum(abs(Pxzp)*dffp),'-m','linewidth',2)
   semilogx(ffp,cumsum(abs(Pyzp)*dffp),'-c','linewidth',2)
   
   plot([.02;10],[0 0],'--k')
   set(gca,'xlim',[.02 10])
   %axis([.02 10 -.02 .02])
   axis('square')
   set(gca,'xtick',[.02 .1 1 10])
   xlabel('Frequency (Hz)')
   ylabel('Cumulative Co-spectra (m^2/s^2)')
   
   subplot(236)
   uncor1 = find(lags <=-timescale);
   uncor2 = find(lags >=timescale);
   uncor = [uncor1(:);uncor2(:)];
   fprintf('number of uncorrelated samples %g\n',length(uncor));
   pct =  100.*length(uncor)/length(xxz);
   fprintf(' xxz is %g percentile.\n',pct)
   if(pct<100),
      hhh=wnormplot_crs(xxz(uncor),1,1,fix(dof));
   end
   hold on
   ax = axis;
   plot([xxz(lag0) xxz(lag0)],ax(3:4),'-r')
   drawnow
   shg
end
rss.pct = pct;
return