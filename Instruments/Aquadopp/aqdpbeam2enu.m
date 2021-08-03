function [u,v,w] = aqdpbeam2enu(ncells,vel_b1,vel_b2,vel_b3,heading,pitch,roll,ornt)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Syntax: 
%   
%  [u,v,w] = aqdpbeam2enu(ncells,vel_b1,vel_b2,vel_b3,heading,pitch,roll,ornt)
%
%  Inputs: 
%          ncells: from hdr file
%          vel_b1: beam 1 velocities
%          vel_b2: beam 2 velocities
%          vel_b3: beam 3 velocities
%          heading: Aquadopp heading (alternatively compass bearing
%                   measured in the field)
%          pitch: Aquadopp pitch
%          roll: Aquadopp roll
%          ornt: from hdr file, or user provided
%
%  Rotate velocities to ENU
%
%  Developed by Dr. Julia C Mullarney, University of Waikato, New Zealand 
%  c. 2013
%
%  Additions made by Benjamin K Norris, 2014
%  Edits: cleaned up and organized original code fragments
%         removed extraneous bits of code
%         13/01/2016: added matlabpool to optimize slow step of code using
%         multicore processing.
%         
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp('Rotating instrument beams in beam to ENU')
pause(1)
T=[6461 -3232 -3232; 0 -5596 5596;1506 1506 1506]/4096; %Transformation matix from hdr file
%now hardwired in to preserve FP numbers

%flip the matrix 
if strcmp(ornt,'DOWN') == 1,
   T(2,:) = -T(2,:);   
   T(3,:) = -T(3,:);    
end
matlabpool(4) %optimize code to use multi-core processing
parfor ii=1:length(pitch)

    hh=pi*(heading-90)/180;
    pp = pi*(pitch(ii))/180;
    rr = pi*(roll(ii))/180;

    % Make heading matrix
    H = [cos(hh) sin(hh) 0; -sin(hh) cos(hh) 0; 0 0 1];

    % Make tilt matrix
    P = [cos(pp) -sin(pp)*sin(rr) -cos(rr)*sin(pp);...
          0             cos(rr)          -sin(rr);  ...
          sin(pp) sin(rr)*cos(pp)  cos(pp)*cos(rr)];

    % Make resulting transformation matrix
    R = H*P*T;

    for nn=1:ncells
        Vbeam = [];
        Vbeam=([[vel_b1(ii,nn)],[vel_b2(ii,nn)],[vel_b3(ii,nn)]]');
        Venu(:,ii,nn)=R*Vbeam;
    end
    
    if (ii/100000)-floor(ii/100000)==0
        disp([num2str(ii) ' samples of ' num2str(numel(pitch)) ' rotated']),
    end
   
end
matlabpool CLOSE

u=squeeze(Venu(1,:,:));
v=squeeze(Venu(2,:,:));
w=squeeze(Venu(3,:,:));

end

