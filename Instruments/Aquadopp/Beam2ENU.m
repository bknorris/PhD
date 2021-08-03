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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Use of this program is described in:
%
% Sullivan, C.M., Warner, J.C., Martini, M.A., Voulgaris, G., 
% Work, P.A., Haas, K.A., and Hanes, D.H. (2006) 
% South Carolina Coastal Erosion Study Data Report for Observations
% October 2003 - April 2004., USGS Open-File Report 2005-1429.
%
% Program written in Matlab v7.1.0 SP3
% Program ran on PC with Windows XP Professional OS.
%
% "Although this program has been used by the USGS, no warranty, 
% expressed or implied, is made by the USGS or the United States 
% Government as to the accuracy and functioning of the program 
% and related program material nor shall the fact of distribution 
% constitute any such warranty, and no responsibility is assumed 
% by the USGS in connection therewith."
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%

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
