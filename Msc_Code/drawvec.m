%
% DRAWVEC draw arrows for vectors
%
% DRAWVEC drawvec(xo,yo,um,vm,atheta,lcolor) is a routine to 
%         draw vectors. This is a fairly low-level routine
%         in that in does no scaling to the vectors.  This function
%         is called primarily by VECPLOT and returns the handle
%         of the vector object drawn.
%
% Inputs: xo,yo  - vector origins; arrow eminates from this point
%         um,vm  - vector magnitudes
%         atheta - arrow-head angle, in degrees
%         lcolor - linecolor , 'r' = red, etc.
%
% Outputs: hp    - the handle to the vector object drawn
%
% Call as:  hp=drawvec(xo,yo,um,vm,atheta,lcolor)
%
% Written by: Brian O. Blanton
%
   function hp=drawvec(xo,yo,um,vm,atheta,lcolor)
   fac = 3.14159/180.;
   arrowtheta = atheta*fac;

% columnate the input vectors to ensure they are 
% column-vectors, not row-vectors
   xo=xo(:);
   yo=yo(:);
   um=um(:);
   vm=vm(:);

% compute and draw arrow shaft
   xe = xo + um;
   ye = yo + vm;
   arrowmag = .25*(sqrt((xo-xe).*(xo-xe)+(yo-ye).*(yo-ye)));
   shafttheta = -atan2((ye-yo),(xe-xo));
   xt = xe-arrowmag.*cos(arrowtheta);
   yt = ye-arrowmag.*sin(arrowtheta);
   x1 = (xt-xe).*cos(shafttheta)+(yt-ye).*sin(shafttheta)+xe;
   y1 = (yt-ye).*cos(shafttheta)-(xt-xe).*sin(shafttheta)+ye;
   xt = xe-arrowmag.*cos(-arrowtheta);
   yt = ye-arrowmag.*sin(-arrowtheta);
   x2 = (xt-xe).*cos(shafttheta)+(yt-ye).*sin(shafttheta)+xe;
   y2 = (yt-ye).*cos(shafttheta)-(xt-xe).*sin(shafttheta)+ye;
   x=ones(length(xo),6);
   y=ones(length(xo),6);
   x=[xo xe x1 xe x2 NaN*ones(size(xo))]';
   y=[yo ye y1 ye y2 NaN*ones(size(yo))]';
   x=x(:);
   y=y(:);
   hp=line(x,y,'LineStyle','-','Color',lcolor);
%
%        Brian O. Blanton
%        Curr. in Marine Sciences
%        15-1A Venable Hall
%        CB# 3300
%        Uni. of North Carolina
%        Chapel Hill, NC
%                 27599-3300
%
%        919-962-4466
%        blanton@marine.unc.edu
%