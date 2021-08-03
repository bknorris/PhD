% vecplot2.m
%
% Usage:
%
% vecplot2(time,uu,vv,date,axis_h)
% time - time vector of julian days
% uu - vector of u (east-west) velocity
% vv - vector of v (north-south) velocity
% time, uu, vv must be same length
% date - 2 element vector of [start_date stop_date] of values to be
%          extracted from time vector
% axis_h current axis handle to be plotted on.

function vecplot2(time,uu,vv,date,axis_h);

timediff = diff(date);
xlim = [(date(1) - timediff*0.02) (date(2) + timediff*0.02)];
scale = [min([min(uu) min(vv)]) max([max(uu) max(vv)])];
vscale = [min(vv) max(vv)];
vscale(1) = min([vscale(1) 0]);
vscale(2) = max([vscale(2) 0]);

%to make several weeks with same axes..
vscale(1)=-25;
vscale(2)=25;

%set(axis_h,'XLim',xlim)
%set(axis_h,'YLim',vscale);
axis([xlim vscale])
old_units = get ( axis_h, 'Units' );
set ( axis_h, 'Units', 'Pixels' );
pos = get ( axis_h, 'Position' );
set ( axis_h, 'Units', old_units );

scale_factor = (diff(xlim)/diff(scale)) * (pos(4)/pos(3));
u_data=scale_factor*uu;
v_data=vv;
x=time;
xp=x;
yp=zeros(size(xp));
xplot=ones(length(xp),2);
yplot=xplot;
xplot(:,1)=x(:);
xplot(:,2)=xp(:)+u_data(:);
xplot(:,3)=x(:);
yplot(:,1)=yp(:);
yplot(:,2)=yp(:)+v_data(:);
yplot(:,3)=yp(:)*nan;
xplot=xplot';
yplot=yplot';
if(~isempty(find(isfinite(u_data(:))==1)))
   graph(1) = plot( date,[0 0] );
   graph(2) = plot ( xplot(:), yplot(:) );...
end
axis([xlim vscale]);
%set(axis_h,'XLim',xlim,'YLim',vscale);
