function  s = uvwstats(u,v,sf,ipost)
% uvwstats - Horizontal statistics of 2-d (e.g. current meter) data
% s = uvwstats(u,v,[[sf],ipost])
%
% u and v are column vectors with x and y data
% sf is an optional length for the axes
% ipost is an optional flag 0=no action, >1=plot input data
 
% Chris Sherwood, USGS
% 9-Jan-2013

if(nargin <4),ipost=0;,end
if(nargin <3),ipost=0;sf=0;,end

mu = mean(u);
mv = mean(v);
C = cov(u,v);
[V,D] = eig(C);

x1 = [.5*sqrt(D(1,1))*V(1,1);-.5*sqrt(D(1,1))*V(1,1)];
y1 = [.5*sqrt(D(1,1))*V(2,1);-.5*sqrt(D(1,1))*V(2,1)];
x2 = [.5*sqrt(D(2,2))*V(1,2);-.5*sqrt(D(2,2))*V(1,2)];
y2 = [.5*sqrt(D(2,2))*V(2,2);-.5*sqrt(D(2,2))*V(2,2)];
[mspd mdir]=pcoord( mu, mv );
[ l1, az1 ] = pcoord( x1(1), y1(1) );
[ l2, az2 ] = pcoord( x2(1), y2(1) );
if(l1 < l2),
  ltemp = l1; aztemp = az1;
  l1 = l2;    az1 = az2;
  l2 = ltemp; az2 = aztemp;
end
sd1 = 2*l1;
sd2 = 2*l2;

s.ubar = mu;
s.vbar = mv;
[S ddir]=pcoord(mu,mv);
s.S = S;
s.az0 = ddir; 
s.sd1 = sd1; % this is also wave std. dev.
s.sd2 = sd2;
s.az1 = az1;
s.az2 = az2;
phid = az1-ddir;
% find acute angle
if( phid < -180 ), phid = phid+180; end
if( phid >  180 ), phid = -(phid-180); end
s.phid = phid;
s.phir = (pi/180)*abs(phid);
s.variables = {...
   'Note: units assume input units are m/s';...
   'Note: angles are geographic (degrees CW from north == v-axis)';...
   's.S = mean speed [m/s]';...
   's.az0 = direction S is toward [degrees]';...
   's.sd1 = std. dev. of major component [m/s]';...
   's.az1 = direction of major component (ambiguous +-180) [degrees]';...
   's.sd2 = std. dev. of minor component [m/s]';...
   's.az2 = direction of minor component (ambiguous +-180) [degrees]';...
   's.phid = az1-az0 (acute angle waves-currents [degrees CW]';...
   's.phir = (pi/180)*abs(phid) [radians, no direction]'};

if(ipost),
   hh1=plot(u(:),v(:),'.b')
   set(hh1,'color',[.9 .4 .4])
   hold on
   axis('square');
   axis([-sf sf -sf sf])
   ts = ['Speed= ',sprintf('%4.2f',mspd),'; Dir.= ',sprintf('%3.0f',mdir) ];
   hh2=text(-.8*sf,.8*sf,ts);
   for i=1:2
      eval(['[ leng(i), az(i) ] = pcoord( x' num2str(i) '(1), y' num2str(i) '(1) );']);
   end;
   ts = ['Major axis: Mag.= ',sprintf('%4.2f',sd1),'; Az.= ',sprintf('%3.0f',az1) ];
   hh2=text(-.8*sf,.6*sf,ts);
   hh3=plot(x1,y1,'-r',x2,y2,'-r','linewidth',2);
   hh4=plot([0; mu],[0; mv],'-b','linewidth',1);
   xlabel('East ({\itu}) component (m/s)')
   ylabel('North ({\itv}) component (m/s)')
   grid
end