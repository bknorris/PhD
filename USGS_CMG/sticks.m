function h = sticks(x,mag,theta,xlims,ylims,magref,u)

% STICKS.M: Creates stick plot.
%			STICKS(X,MAG,THETA) makes a stick plot in the current axis where
%			X is the position along the x-axis where each stick starts,
%			MAG is the magnitude and THETA is the compass direction (0 is in
%			the direction of the y-axis).
%			For multiple plots in the same axis, X, MAG and THETA are MxN
%			matrices where each of the N columns represent one plot. Should
%			any plot have less values than the plot with the most values, the
%			superfluous matrix elements can be filled with NaNs. The plots are
%			from the bttom of the axis and upwards.
%			Optional arguments: H = STICKS(X,MAG,THETA,XLIMS,YLIMS,MAGREF,UNITS)
%			XLIMS sets the limits of the x-axis (default: [min(X) max(X)]).
%			YLIMS does the same for y-axis, but note that the y-axis limits
%			scale the stick vector magnitudes (default: [-max(mag) max(mag)]).
%			For multiple plots, YLIMS is multiplied by N.
%			MAGREF is the length of an optional reference vector. Default is no
%			reference vector.
%			UNIT can be set to 'rad' if the direction is given in radians
%			instead of degrees (the default).
%			Default values for optional arguments are indicated by the
%			empty matrix ([]).
%			H contains the handles to all the sticks (line objects). If MAGREF
%			has been set the last handle in H refers to the reference vector.
%
%			NOTE! It is presently not possible to change the position or size
%			of the axes once the plot has been drawn, as this will distort the
%			size and direction of the sticks. Instead, delete the axes object,
%			create a new at the correct position and then run STICKS again.

%			Olof Liungman, 1999
%			Dept. of Oceanography, Earth Sciences Centre
%			Göteborg University, Sweden
%			E-mail: olof.liungman@oce.gu.se

if nargin<3
  error('Not enough input arguments!')
elseif nargin==3
  xlims = [];
  ylims = [];
  magref = [];
  u = 'deg';
elseif nargin==4
  ylims = [];
  magref = [];
  u = 'deg';
elseif nargin==5
  magref = [];
  u = 'deg';
elseif nargin==6
  u = 'deg';
end

if ~strcmp(u,'rad')&~strcmp(u,'deg')
  error('Units input argument must be ''deg'' or ''rad''!')
end

sx = size(x); smag = size(mag); stheta = size(theta);
if (min(sx)+min(smag)+min(stheta))==3
	if (max(sx)~=max(smag))|(max(smag)~=max(stheta))
		error('Input vectors must be of equal length!')
	end
	x = x(:); mag = mag(:); theta = theta(:);
  [m,nplots] = size(x);
else
	if (any(sx~=smag))|(any(smag~=stheta))
  	error('Input matrices must be of equal size!')
	end
	[m,nplots] = size(x);
end

if strcmp(u,'deg')
  theta = pi*theta/180;
end

if isempty(xlims)
  xlims = [min(min(x)) max(max(x))];
end
if isempty(ylims)
  ylims = [-1 1]*max(max(mag));
end
ylims = nplots*ylims;

ah = gca;

set(ah,'YLim',ylims,'XLim',xlims,'Units','pixels',...
		'Box','on')

axes_pos = get(ah,'Position');
pixels_per_yunit = axes_pos(4)/(ylims(2)-ylims(1));
pixels_per_xunit = axes_pos(3)/(xlims(2)-xlims(1));
xdata_ratio = pixels_per_yunit/pixels_per_xunit;

handles = [];
for n = 1:nplots

  xn = [x(:,n) mag(:,n).*sin(theta(:,n))*xdata_ratio+x(:,n)];
	yoffset = ylims(1)+(n-0.5)/nplots*diff(ylims);
  yn = [zeros(m,1) mag(:,n).*cos(theta(:,n))]+yoffset;
	handlesn = plotsticks(xn,yn,[min(x(:,n)) max(x(:,n))],...
						 						[1 1]*yoffset,ah);
	handles = [handles;handlesn];

end

if ~isempty(magref)
	refvec_pos = [xlims(1)+diff(xlims)/10 ylims(1)+diff(ylims)/10];
	refhandle = line(refvec_pos(1)+[0 magref]*xdata_ratio,...
									 refvec_pos(2)+[0 0],'Color','b','LineStyle','-');
	handles = [handles;refhandle];
end

if nargout==1
  h = handles;
end

%----------------------------------------------------------------
function h = plotsticks(x,y,x0,y0,axeshandle)

h = line(x',y','Color','b','LineStyle','-');
xaxescolor = get(axeshandle,'XColor');
line(x0,y0,'Color',xaxescolor,'LineStyle','-');