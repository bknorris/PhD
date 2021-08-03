function out = cmg_pca(u,v,s,ipost)
% PCA  Principal components of 2-d (e.g. current meter) data
%
% function pca(u,v,[[s],ipost])
%
% u and v are column vectors with x and y data
% s is an optional length for the axes
% ipost is an optional flag 0=no action, 1=plot input data
% out is a structure containing the fields:
% mspd = mean speed
% mdir = mean direction
% ax1 = [axis1 length & azimuth]
% ax2 = [axis2 length & azimuth]
 
% Chris Sherwood, USGS
% March 17, 1999

if(nargin <4),ipost=0;end
if(nargin <3),ipost=0;s=0;end

mu = mean(u);
mv = mean(v);
C = cov(u,v);
[V,D] = eig(C);

x1 = [.5*sqrt(D(1,1))*V(1,1);-.5*sqrt(D(1,1))*V(1,1)];
y1 = [.5*sqrt(D(1,1))*V(2,1);-.5*sqrt(D(1,1))*V(2,1)];
x2 = [.5*sqrt(D(2,2))*V(1,2);-.5*sqrt(D(2,2))*V(1,2)];
y2 = [.5*sqrt(D(2,2))*V(2,2);-.5*sqrt(D(2,2))*V(2,2)];

% clf
% plot(x1,y1,'-r',x2,y2,'-r')
% axis('square')
% if(s==0),s=ceil(10*max([max([x1,y1,x2,y2]),mu,mv]))/10;end
% if s==1,axis([-s s -s s]);end
% xlabel('U component')
% ylabel('V component')
% grid
% hold on
if(ipost==1),
  plot(u,v,'.r')
end
[mspd,mdir]=pcoord( mu, mv );
% plot([0; mu],[0; mv],'-b')
% hold on
% fprintf( 'Mean u = %f    Mean v = %f\n', mu, mv );
% fprintf( 'Mean magnitude: %f   Direction: %f\n', mspd,mdir );
[ l1, az1 ] = pcoord( x1(1), y1(1) );
[ l2, az2 ] = pcoord( x2(1), y2(1) );
% fprintf( 'Axis 1 Length: %f   Azimuth: %f\n', l1*2, az1 );
% fprintf( 'Axis 2 Length: %f   Azimuth: %f\n', l2*2, az2 );
out.mspd = mspd;out.mdir = mdir;out.ax1 = [l1 az1];out.ax2 = [l2 az2];


