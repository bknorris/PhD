


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Interpolate (with spline) met station data to Aquadopp times and subtract air pressure
%  variations.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

 
load air_pressure_raw
airpressure=spline(yd,airpressure_dbar,aquadopp.yearday);
if aquadopp.yearday(end)>=yd(end)  %cut to same time
   [junk indap]=min(abs(yd(end)-aquadopp.yearday));
   airpressure(indap:length(aquadopp.pressure))=airpressure(indap);
end
if aquadopp.yearday(1)<=yd(1)
   [junk indap]=min(abs(yd(1)-aquadopp.yearday));
   airpressure(1:indap)=airpressure(indap);
end


offset=[0 0 0 0 0.35 0.35 0.35 0.37];
aquadopp.pressure_corrected=aquadopp.pressure+mean(airpressure)-airpressure;  %-offset(adnumber)
%%%BEN -- CHECK THIS IS THE RIGHT WAY ROUND
%Also - I've removed the offset from this bit as it's better done after the
%temperature correction anyway

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  End air pressure correction
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% T-P varation
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure(1)
    plot(aquadopp.yearday,aquadopp.pressure_corrected)
    hold on


nsec=2;  %Number of sections... probably 2 will be right (before inundation and after indunatio)

for ii=1:nsec  %click start and end bits of 'out of water sections'
    
    clear q
    q=ginput(2);
    
    [junk startout(ii)]=min(abs(q(1,1)-aquadopp.yearday));
    [junk endout(ii)]=min(abs(q(2,1)-aquadopp.yearday));
    
end %endfor
    
figure(2)
hold on


for ii=1:nsec   %Combines all data into single variable to plot it.
    if ii==1
        p_all=aquadopp.pressure_corrected(startout(ii):endout(ii))';
        temp_all=aquadopp.temperature(startout(ii):endout(ii))';
    else
        p_all=[p_all pressure_corrected(startout(ii):endout(ii))'];
        temp_all=[temp_all aquadopp.temperature(startout(ii):endout(ii))'];
    end
    plot(aquadopp.temperature(startout(ii):endout(ii)),pressure_corrected(startout(ii):endout(ii)),'o')
end

p=polyfit(temp_all,p_all,1); %Fit and plot linear relationship 
xlim=get(gca,'xlim');
plot(xlim,p(1)*xlim+p(2),'k')


new_corrected_pressure=pressure_corrected-p(1)*aquadopp.temperature; %Removes T-variation
figure(1)
    plot(aquadopp.year,new_corrected_pressure,'m')  %  See how this sits... will give you a better idea of offset...
    
offset=[0 0 0 0 0.35 0.35 0.35 0.37];

%%%%CEHCK THIS AND RENAME TO SOMETHING SENSIBLE!!
% Suggest ending up with final variable aquadopp.pressure_corrected and
% call the interim bits something else.












%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Rotate velocities to xyz instrument coordinates
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Tmat=[6461 -3232 -3232; 0 -5596 5596;1506 1506 1506]/4096; %Transformation matix from hdr file

%if statusbit0 == 1, downward looking
 %  Tmat(2,:) = -Tmat(2,:);   
 %  Tmat(3,:) = -Tmat(3,:);   
%end

for ii=1:ncells
    clear beamvels
    beamvels=([[aquadopp.vel_b1(:,ii)],[aquadopp.vel_b2(:,ii)],[aquadopp.vel_b3(:,ii)]]');
    xyz_vels(:,:,ii)=Tmat*beamvels;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Transform velocities to ENU coordinates
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


for ii=1:length(aquadopp.heading)

    clear hh pp rr H P R
    
     
    hh=pi*(VALUE_FROM_NOTEBOOK-90)/180;  %%%%*****  BEN - here value is the value we measured - or you can use aquadopp.heading(ii) if it's OK
    
    pp = pi*(aquadopp.pitch(ii))/180;
    rr = pi*(aquadopp.roll(ii))/180;

    % Make heading matrix
    H = [cos(hh) sin(hh) 0; -sin(hh) cos(hh) 0; 0 0 1];

    % Make tilt matrix
    P = [cos(pp) -sin(pp)*sin(rr) -cos(rr)*sin(pp);...
          0             cos(rr)          -sin(rr);  ...
          sin(pp) sin(rr)*cos(pp)  cos(pp)*cos(rr)];

    % Make resulting transformation matrix
    R = H*P*Tmat;

    for nn=1:ncells
        clear beamvels
        beamvels=([[aquadopp.vel_b1(ii,nn)],[aquadopp.vel_b2(ii,nn)],[aquadopp.vel_b3(ii,nn)]]');
        ENU_vels(:,ii,nn)=R*beamvels;
    end

    if (ib/5000)-floor(ii/5000)==0  % monitor progress of rotation to stop Julia getting impatient
       disp('ii '),ii 
    end
   
end


u=squeeze(ENU_vels(1,:,:));
v=squeeze(ENU_vels(2,:,:));
w=squeeze(ENU_vels(3,:,:));


%%%  Then we need to do magnetic correction....


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Rotate velocitites so they are relative to true N not magnetic North
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% Rotate to True North coordinates
mag_offset =VALUE*pi/180; % degrees   %VALUE=??

u_old=u;
v_old=v;


u= u_old.*(ones(size(u_old))*cos(mag_offset)) + ...
         v_old.*(ones(size(v_old))*sin(mag_offset));
v=-u_old.*(ones(size(u_old))*sin(mag_offset)) + ...
       v_old.*(ones(size(v_old))*cos(mag_offset));

clear u_old v_old
