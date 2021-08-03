function Venu=Beam2ENU(Vbeam,heading,pitch,roll,ornt);
% function Venu=Beam2ENU(Vbeam,T,heading,pitch,roll,up)
%
% Beam2ENU transforms beam data to ENU coordinates. 
% This version assumes that the heading, pitch and roll are fixed.
%
% Vbeam is the 3D velocity in beam coordinates. Array can be in 
% rows or columns as long as the number of data points is >3. 
% 
% T is the transformation matrix for beam to xyz coordinates, 
% Examples 
% T = [6461 -3232 -3232; 0 -5596 5596; 1506 1506 1506]; % profiler 
% T = [2896 2896 0; -2896 2896 0; -2896 -2896 5792]; % std current meter 
% 
% heading, pitch and roll are all in degrees 
% ornt==0 if the instrument points up, ornt==1 if down 

%
T =[6461       -3232       -3232 ;...
      0       -5596        5596 ;....
    1506        1506        1506];

%
T = T/4096;   % Scale the transformation matrix correctly to floating point numbers
%
if ornt == 1,
   T(2,:) = -T(2,:);   
   T(3,:) = -T(3,:);   
end
%
% Put data in columns for the calculation
[nr,nc] = size( Vbeam );                 % find size of x
car=1;
if(nr<nc),
  Vbeam=Vbeam';
  car=0;
  end;


Vxyz = T*Vbeam;
%hh = pi*heading/180;
hh = pi*(heading-90)/180;                      %comment here
%  The 90 deg offset is because Vx points north and Vy points 
%  west when heading==0.
pp = pi*pitch/180;
rr = pi*roll/180;

% Make heading matrix
H = [cos(hh) sin(hh) 0; -sin(hh) cos(hh) 0; 0 0 1];

% Make tilt matrix
P = [cos(pp) -sin(pp)*sin(rr) -cos(rr)*sin(pp);...
      0             cos(rr)          -sin(rr);  ...
      sin(pp) sin(rr)*cos(pp)  cos(pp)*cos(rr)];

% Make resulting transformation matrix
R = H*P;  

Venu = R*Vxyz;

% return array to column/row format of Vbeam
if (car==0), Venu=Venu'; end;

