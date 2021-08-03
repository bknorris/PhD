function [x, y] = sticksafe(x0, y0, theLength, theAngle)

% sticksafe -- Intelligent, oriented sticks.
%  sticksafe('demo') demonstrates itself with sticks centered
%   on the origin, each one-unit long.
%  sticksafe(N) demonstrates itself with N sticks.
%  sticksafe(x0, y0, theLength, theAngle) draws sticks that
%   start at (x0, y0), with theLength (in y-axis units), and
%   theAngle (degrees, counter-clockwise from +x).  The
%   variables should be the same size, but any can be a scalar,
%   just so long as the x0 and/or y0 array is the full size.
%   The "ResizeFcn" of the figure is set to update the sticks
%   automatically.  The "Tag" of each stick is the mfilename.
%   Properties of the arrows, such as "color", can be adjusted
%   at any time.  The arrows are Matlab "patch" objects.
%  h = sticksafe(...) draws the sticks and returns the handle.
%  [x, y] = sticksafe(...) returns the (x, y) data for such
%   sticks, one row per arrow, but does not draw them.
%  sticksafe (no argument) redraws existing sticks. This is
%   useful whenever the window is resized or the x or y limits
%   change.  (Recall that printing causes the "ResizeFcn" to
%   be called twice.)
%
% Note: this routine leaves the axes in "manual" mode.  Use
%  "axis auto" to revert to automatic axis limits.
 
% Copyright (C) 2000 Dr. Charles R. Denham, ZYDECO.
%  All Rights Reserved.
%   Disclosure without explicit written consent from the
%    copyright owner does not constitute publication.
 
% Version of 12-Jan-2000 14:22:59.
% Updated    27-Jul-2001 15:01:56.

RCF = 180 / pi;

% Resize.

if nargin < 1
	oldGCA = gca;
	h = findobj(gcf, 'Type', 'line', 'Tag', mfilename);
	for i = 1:length(h)
		p = get(h(i), 'Parent');
		axes(p)
		u = get(h(i), 'UserData');
		[xx, yy] = feval(mfilename, u(:, 1), u(:, 2), u(:, 3), u(:, 4));
		set(h(i), 'XData', xx(:), 'YData', yy(:));
	end
	axes(oldGCA)
	return
end

% Demonstration.

if nargin == 1
	if isequal(x0, 'demo')
		help(mfilename)
		x0 = 16;
	elseif ischar(x0)
		x0 = eval(x0);
	end
	theName = [mfilename ' demo'];
	f = findobj('Type', 'figure', 'Name', theName);
	if ~any(f)
		f = figure('Name', theName);
	end
	figure(max(f))
	delete(get(f, 'Children'))
	n = max(1, round(x0));
	x0 = zeros(1, n);
	ang = linspace(0, 360, length(x0)+1);
	ang(end) = [];
	h = feval(mfilename, x0, 0, 1, ang);
	set(gca, 'Xlim', [-n n], 'YLim', [-2 2])
	feval(mfilename)
	set(gcf, 'WindowButtonDownFcn', ...
			['if zoomsafe(''down''), ' mfilename ', end'])
	if nargout > 0, x = h; end
	return
end

% Initialize.

if nargin > 1
	if nargout == 2, x = []; y = []; end
	if length(x0) == 1
		x0 = x0 * ones(size(y0));
	elseif length(y0) == 1
		y0 = y0 * ones(size(x0));
	end
	x0 = reshape(x0, [1 prod(size(x0))]);
	y0 = reshape(y0, size(x0));
	if nargin < 3, theLength = 1; end
	if nargin < 4, theAngle = 0; end
	if length(theLength) == 1
		theLength = theLength * ones(size(x0));
	end
	if length(theAngle) == 1
		theAngle = theAngle * ones(size(x0));
	end

	theLength = reshape(theLength, size(x0));
	theAngle = reshape(theAngle, size(x0));

	axes(gca)
	oldUnits = get(gca, 'Units');
	set(gca, 'Units', 'pixels')
	thePosition = get(gca, 'Position');
	set(gca, 'Units', oldUnits)
	theWidth = thePosition(3);   % pixels.
	theHeight = thePosition(4);   % pixels.
	
	axis('manual')
	dx = diff(get(gca, 'XLim'));
	dy = diff(get(gca, 'YLim'));
	dydx = dy / dx;   % Not used.
	dxdp = dx / theWidth;   % sci/pixel.
	dydp = dy / theHeight;   % sci/pixel.
	scale = dxdp / dydp;

	zz = exp(sqrt(-1) .* theAngle ./ RCF) .* theLength;
	
	xx = zeros(3, prod(size(x0)));
	xx(1, :) = x0;
	xx(2, :) = x0 + real(zz)*scale;
	xx(3, :) = NaN;
	
	yy = zeros(size(xx));
	yy(1, :) = y0;
	yy(2, :) = y0 + imag(zz);
	yy(3, :) = NaN;
	
	if nargout < 2
		h = line(xx(:), yy(:), ...
			'Tag', mfilename, ...
			'UserData', [x0(:) y0(:) theLength(:) theAngle(:)]);
		if nargout > 1, x = h; end
	elseif nargout == 2
		x = xx;
		y = yy;
	end
	set(gcf, 'ResizeFcn', mfilename)
	
	if nargout == 1, x = h; end
end