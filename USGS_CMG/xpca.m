function [s]=xpca(u,v,ipost)
% Principal components of 2-d (e.g. current meter) data
%
% function pca(u,v,[[s],ipost])
%
% u and v are column vectors with x and y data
% s is an optional length for the axes
% ipost is an optional flag 0=scatter only, 1=ellips only, 2=both
 
% Chris Sherwood, USGS
% March 17, 1999

% Modified by jpx on 6-14-99 to be used in CMGTooL
% Modified by jpx on 10-14-99 to plot ellipses

mu = mean(u);
mv = mean(v);
C = cov(u,v);
[V,D] = eig(C);

x1 = [.5*sqrt(D(1,1))*V(1,1);-.5*sqrt(D(1,1))*V(1,1)];
y1 = [.5*sqrt(D(1,1))*V(2,1);-.5*sqrt(D(1,1))*V(2,1)];
x2 = [.5*sqrt(D(2,2))*V(1,2);-.5*sqrt(D(2,2))*V(1,2)];
y2 = [.5*sqrt(D(2,2))*V(2,2);-.5*sqrt(D(2,2))*V(2,2)];
[mspd mdir]=xpcoord( mu, mv );


% clf
oldunits=get(gcf,'units');

set(gcf,'units','norm');

hh2(1)=subplot(212);
set(hh2(1),'position',[.1 .05 .8 .19],'visible','off');

hh2(2)=text(.2,.8,['Speed= ' num2str(mspd) ';  Direction= ' num2str(mdir) ] );
set(hh2(2),'color','b','visible','off');

for i=1:2
	eval(['[ leng(i), az(i) ] = xpcoord( x' num2str(i) '(1), y' num2str(i) '(1) );']);
end;

if leng(1)>=leng(2)
	majr=leng(1);majaz=az(1);
	minr=leng(2);minaz=az(2);
else
	majr=leng(2);majaz=az(2);
	minr=leng(1);minaz=az(1);
end;

hh2(3)=text(.2,.5,['Major axis: Magnitude= ' num2str(majr*2) ';  Azimuth= ' num2str(majaz) ] );
set(hh2(3),'color','r','visible','off');
hh2(4)=text(.2,.2,['Minor axis: Magnitude= ' num2str(minr*2) ';  Azimuth= ' num2str(minaz) ] );
set(hh2(4),'color','r','visible','off');
if ipost>0
	set(hh2,'visible','on');
axis off;
end;

hand=subplot(211);
if ipost~=1
	ticks=ceil(length(u)/5000);
% 	plot(u(1:ticks:end),v(1:ticks:end),'.g','markersize',1)
	xx=round(length(u)/ticks);
	indx=round(rand(xx,1)*length(u));
	indx(indx==0)=[];
	plot(u(indx),v(indx),'.g','markersize',1);
end;
hold on;
set(hand,'position',[.1 .35 .8 .6]);

hh1=plot(x1,y1,'-r',x2,y2,'-r','linewidth',2);
% set(hh(1),'linewidth',3);
% set(hh(2),'linewidth',3);

theta=2*pi*(1:64)/64;
theta=[0 theta];
xx=majr*cos(theta);
yy=minr*sin(theta);
angle=-majaz*pi/180+pi/2;
xxx=xx*cos(angle)-yy*sin(angle);
yyy=xx*sin(angle)+yy*cos(angle);
hh1(3)=plot(xxx,yyy,'r');


axis('square');

%if ipost>0
	s=min([max(abs(u)),max(abs(v))]);
	s=ceil(0.75*s/10)*10;
%else
%	s=ceil(max([max([x1,y1,x2,y2]),mu,mv])/10)*10;
%end;
axis([-s s -s s])
% xlabel('U component')
% ylabel('V component')
% grid
hh1(4)=plot([0; mu],[0; mv],'-b','linewidth',1);
if ipost<1
	set(hh1,'visible','off');
end;
box on;
set(gcf,'units',oldunits);