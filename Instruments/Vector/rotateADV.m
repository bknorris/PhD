function ADV = rotateADV(ADV,metadata)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   A function that transforms Nortek Vector velocities from either BEAM or 
%   XYZ to ENU coordinates. This function also applies an adjustment for 
%   the supplied magnetic declination.
%
%   ADV - a structure ADV containing velocities and sensor data
%   Metadata - user defined metadata. For a  complete list of required 
%   metadata fields, see adv_dataprocess.m
%
%   Dependencies: despikeADV.m
%
%   Contains a partial adaptation of vec2ncBurstNative.m
%   Written by Kurt J Rosenberger
%   krosenberger@usgs.gov
%   USGS Pacific Coastal Marine Science Center
%
%   This script was compiled by Benjamin K Norris, 2015
%   University of Waikato, New Zealand
%  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('Performing transformations to Earth Coordinates...')
T = ADV.Metadata.instmeta.transfm;
samprate = ADV.Metadata.instmeta.samprate;
% NOTE: For the Vector the instrument is defined to be in 
%       the UP orientation when the probe is pointing down.
%for CABLED vectors only, the z axis points in the same direction as the
%case
switch metadata.orientation
    case 'UP' %probe points UP, z-axis points DOWN
        ADV.Metadata.velocity_height = ADV.Metadata.probe_height + 150;
        %in this case, we need to change the sign of rows 2 and 3 in the
        %transformation matrix
        T(2,:) = -T(2,:);
        T(3,:) = -T(3,:);
        if strcmp(ADV.Metadata.instmeta.coordsys,'XYZ')
            ADV.V2 = ADV.V2*-1;ADV.V3 = ADV.V3*-1;
            fprintf('Coordinate system is XYZ and the user states that the probe is up\nChanging sign of Y and Z directions\n')
        end
    case 'HORIZONTAL'
        if strcmp(metadata.vector_type,'FIXED')
            ADV.Metadata.velocity_height = ADV.Metadata.probe_height;
        elseif strcmp(metadata.vector_type,'FLEXIBLE')
            if strcmp(metadata.probe_orientation,'UP') && strcmp(metadata.case_rotation,'UP') %case and probe are up, z-axis points up
                ADV.Metadata.velocity_height = ADV.Metadata.probe_height + 150;
            elseif strcmp(metadata.probe_orientation,'HORIZONTAL') && strcmp(metadata.case_rotation,'HORIZONTAL')
                ADV.Metadata.velocity_height = ADV.Metadata.probe_height;
                warning('User says that the probe and case are horizontal. This may cause issues with rotations')
            elseif strcmp(metadata.probe_orientation,'DOWN') && strcmp(metadata.case_rotation,'DOWN') %case and probe point down, z-axis points down, flip Y,Z and T
                ADV.Metadata.velocity_height = ADV.Metadata.probe_height - 150;
                %in this case, we need to change the sign of rows 2 and 3 in the
                %transformation matrix
                T(2,:) = -T(2,:);
                T(3,:) = -T(3,:); 
                if strcmp(ADV.Metadata.instmeta.coordsys,'XYZ')
                    ADV.V2 = ADV.V2*-1;ADV.V3 = ADV.V3*-1;
                    fprintf('Coordinate system is XYZ and the user states that the probe and case are down\nChanging sign of Y and Z directions\n')
                end
            elseif strcmp(metadata.probe_orientation,'DOWN') && strcmp(metadata.case_rotation,'UP') %case is up probe is down, z-axis is up
                ADV.Metadata.velocity_height = ADV.Metadata.probe_height - 150;
            elseif strcmp(metadata.probe_orientation,'UP') && strcmp(metadata.case_rotation,'DOWN') %case is down but probe is up, flip Y,Z and T
                ADV.Metadata.velocity_height = ADV.Metadata.probe_height + 150;
                %in this case, we need to change the sign of rows 2 and 3 in the
                %transformation matrix
                T(2,:) = -T(2,:);
                T(3,:) = -T(3,:);
                if strcmp(ADV.Metadata.instmeta.coordsys,'XYZ')
                    ADV.V2 = ADV.V2*-1;ADV.V3 = ADV.V3*-1;
                    fprintf('Coordinate system is XYZ and the user states that the probe is up and the case is down\nChanging sign of Y and Z directions\n')
                end
            end
        end
        if isfield(metadata,'b1')
            msg = 'Beams have been swapped for rotations';
            ADV.Metadata.swap_beams = msg;fprintf([msg '\n'])
            metadata.swap_beams = msg;
        end
    case 'DOWN' %probe points DOWN, z-axis points UP
        ADV.Metadata.velocity_height = ADV.Metadata.probe_height - 150;
        fprintf('Coordinate system is XYZ and the user states that the probe is down')
end

%correct heading for magnetic variation, throw error if this value is not
%supplied
try
    magvardeg = metadata.magdec;
catch
end
if isempty(magvardeg)
    disp('WARNING: no magnetic declination value is supplied')
    metadata.magdec = input('Enter magnectic declination value to be applied: ');
end

%preallocate velocity variables, define heading pitch and roll
[M,N] = size(ADV.V1);
ADV.U = zeros(M,N);ADV.V = zeros(M,N);ADV.W = zeros(M,N);
if isfield(metadata,'fixed_heading') %user defined fixed heading value
    disp(['User requests a fixed heading of ' num2str(metadata.fixed_heading) ' degrees']);
    heading = ones(M,1)*metadata.fixed_heading;
else
    heading = interp(ADV.Sensor.Heading,samprate);
end
if metadata.zero_pitch %if user wants zero'd pitch/roll
    pitch = zeros(M,1);
else
    pitch = interp(ADV.Sensor.Pitch,samprate);
end
if metadata.zero_roll
    roll = zeros(M,1);
else
    roll = interp(ADV.Sensor.Roll,samprate);
end
V1 = ADV.V1;V2 = ADV.V2;V3 = ADV.V3;

%adjust heading by the magdec supplied
fprintf('Heading is being rotated by %f\n',magvardeg);
heading = heading + magvardeg;
heading(heading >= 360) = heading(heading >= 360)-360;
heading(heading < 0) = heading(heading < 0)+360;
magvar = magvardeg(1)*pi/180;

%we must differentiate between vectors with fixed and flexible heads. fixed
%head vectors will have velocity measurements that are consistent with
%heading/pitch/roll. flexible head vectors may have head orientations that
%are different from the case, and therefore need special circumstances
%for every possible configuration of head and case. 

%define beam velocities
switch metadata.vector_type
    case 'FIXED'
        if isfield(metadata,'swap_beams') %if user requests to swap beam velocities for transformations
            b1 = metadata.b1;b2 = metadata.b2;b3 = metadata.b3;
            if strcmp(b1,'-b2')
                V1 = ADV.V2.*-1;
            end
            if strcmp(b1,'-b3')
                V1 = ADV.V3.*-1;
            end
            if strcmp(b1,'b2')
                V1 = ADV.V2;
            end
            if strcmp(b1,'b3')
                V1 = ADV.V3;
            end
            if strcmp(b2,'-b1')
                V2 = ADV.V1.*-1;
            end
            if strcmp(b2,'-b3')
                V2 = ADV.V3.*-1;
            end
            if strcmp(b2,'b1')
                V2 = ADV.V1;
            end
            if strcmp(b2,'b3')
                V2 = ADV.V3;
            end
            if strcmp(b3,'-b1')
                V3 = ADV.V1.*-1;
            end
            if strcmp(b3,'-b2')
                V3 = ADV.V2.*-1;
            end
            if strcmp(b3,'b1')
                V3 = ADV.V1;
            end
            if strcmp(b3,'b2')
                V3 = ADV.V2;
            end
        end
        switch ADV.Metadata.instmeta.coordsys
            case 'ENU'
                disp('Data are already in Earth Coordinates, applying magnetic variation')
                ADV.U =  V1*cos(magvar) + V2*sin(magvar);
                ADV.V = -V1*sin(magvar) + V2*cos(magvar);
                ADV.W = V3;
            case 'XYZ'
                disp('Data are in XYZ coordinates - applying rotations')
                for i = 1:M
                    hh = pi*(heading(i)-90)/180;
                    pp = pi*pitch(i)/180;
                    rr = pi*roll(i)/180;
                    
                    H = [cos(hh) sin(hh) 0; -sin(hh) cos(hh) 0; 0 0 1];
                    
                    % Make tilt matrix
                    P = [cos(pp) -sin(pp)*sin(rr) -cos(rr)*sin(pp);...
                        0             cos(rr)          -sin(rr);  ...
                        sin(pp) sin(rr)*cos(pp)  cos(pp)*cos(rr)];
                    % Make resulting transformation matrix depending on coord system
                    R = H*P;
                    %transform the data using rotated heading to correct for magnetic variation
                    % do each sample individually for now
                    V = [V1(i,:) V2(i,:) V3(i,:)]';
                    vel = R*V;
                    ADV.U(i,:) = vel(1);
                    ADV.V(i,:) = vel(2);
                    ADV.W(i,:) = vel(3);
                if (i/100000)-floor(i/100000)==0
                    disp([num2str(i) ' samples of ' num2str(numel(heading)) ' rotated']),
                end
                end
            case 'BEAM'
                disp('Data are in BEAM coordinates - applying beam transformations')
                for i = 1:M
                    hh = pi*(heading(i)-90)/180;
                    pp = pi*pitch(i)/180;
                    rr = pi*roll(i)/180;
                    
                    H = [cos(hh) sin(hh) 0; -sin(hh) cos(hh) 0; 0 0 1];
                    
                    % Make tilt matrix
                    P = [cos(pp) -sin(pp)*sin(rr) -cos(rr)*sin(pp);...
                        0             cos(rr)          -sin(rr);  ...
                        sin(pp) sin(rr)*cos(pp)  cos(pp)*cos(rr)];
                    % Make resulting transformation matrix depending on coord system
                    R = H*P*T;
                    %transform the data using rotated heading to correct for magnetic variation
                    % do each sample individually for now
                    V = [V1(i,:) V2(i,:) V3(i,:)]';
                    vel = R*V;
                    ADV.U(i,:) = vel(1);
                    ADV.V(i,:) = vel(2);
                    ADV.W(i,:) = vel(3);
                    if (i/100000)-floor(i/100000)==0
                        disp([num2str(i) ' samples of ' num2str(numel(heading)) ' transformed']),
                    end
                end
        end
    case 'FLEXIBLE'
        %Contingencies:
        %CASE1: if the case and the sensor are pointing in the same direction:
        %rotate using the heading provided, and pitch/roll. 
        %CASE2: if the case and sensor head do not point in the same direction:
        %apply a rotation to the heading such that probe_heading = case_heading.
        if isfield(metadata,'swap_beams') %if user requests to swap beam velocities for transformations
            b1 = metadata.b1;b2 = metadata.b2;b3 = metadata.b3;
            if strcmp(b1,'-b2')
                V1 = ADV.V2.*-1;
            end
            if strcmp(b1,'-b3')
                V1 = ADV.V3.*-1;
            end
            if strcmp(b1,'b2')
                V1 = ADV.V2;
            end
            if strcmp(b1,'b3')
                V1 = ADV.V3;
            end
            if strcmp(b2,'-b1')
                V2 = ADV.V1.*-1;
            end
            if strcmp(b2,'-b3')
                V2 = ADV.V3.*-1;
            end
            if strcmp(b2,'b1')
                V2 = ADV.V1;
            end
            if strcmp(b2,'b3')
                V2 = ADV.V3;
            end
            if strcmp(b3,'-b1')
                V3 = ADV.V1.*-1;
            end
            if strcmp(b3,'-b2')
                V3 = ADV.V2.*-1;
            end
            if strcmp(b3,'b1')
                V3 = ADV.V1;
            end
            if strcmp(b3,'b2')
                V3 = ADV.V2;
            end
        end
        if strcmp(metadata.vector_type,'FLEXIBLE')
            if ~strcmp(metadata.orientation,'HORIZONTAL')
                fprintf('WARNING: current version of rotateADV.m does not support non-horizontal orientations\nof flexible vectors\n')
                return
            end
        end
        if ~isfield(metadata,'probe_heading')
            disp('User claims the instrument is type FLEXIBLE - missing probe heading information')
            metadata.probe_heading = input('Enter probe heading: ');
        end
        hr = range(heading)/2;hm = mean(heading);h1 = floor(hm-hr);h2 = ceil(hm+hr);
        if ~(metadata.probe_heading >= h1 && metadata.probe_heading <= h2)
            probe_heading = ones(M,1)*metadata.probe_heading;
            heading = probe_heading+magvardeg; %keep the same range as the original heading
            heading(heading >= 360) = heading(heading >= 360)-360;
            heading(heading < 0) = heading(heading < 0)+360;
        end
        if ~isfield(metadata,'case_rotation')
            disp('User claims the instrument is type FLEXIBLE - missing case rotation information')
            metadata.case_rotation = input('Enter case rotation [UP/DOWN]: ');
        end       
        switch ADV.Metadata.instmeta.coordsys
            case 'ENU'
                disp('Data are already in Earth Coordinates, applying magnetic variation')
                ADV.U =  V1*cos(magvar) + V2*sin(magvar);
                ADV.V = -V1*sin(magvar) + V2*cos(magvar);
                ADV.W = V3;
            case 'XYZ'
                disp('Data are in XYZ coordinates - applying rotations')
                for i = 1:M
                    hh = pi*(heading(i)-90)/180;
                    pp = pi*pitch(i)/180;
                    rr = pi*roll(i)/180;
                    
                    H = [cos(hh) sin(hh) 0; -sin(hh) cos(hh) 0; 0 0 1];
                    
                    % Make tilt matrix
                    P = [cos(pp) -sin(pp)*sin(rr) -cos(rr)*sin(pp);...
                        0             cos(rr)          -sin(rr);  ...
                        sin(pp) sin(rr)*cos(pp)  cos(pp)*cos(rr)];
                    % Make resulting transformation matrix depending on coord system
                    R = H*P;
                    %transform the data using rotated heading to correct for magnetic variation
                    % do each sample individually for now
                    V = [V1(i,:) V2(i,:) V3(i,:)]';
                    vel = R*V;
                    ADV.U(i,:) = vel(1);
                    ADV.V(i,:) = vel(2);
                    ADV.W(i,:) = vel(3);
                    if (i/100000)-floor(i/100000)==0
                        disp([num2str(i) ' samples of ' num2str(numel(heading)) ' rotated']),
                    end
                end
            case 'BEAM'
                disp('Data are in BEAM coordinates - applying beam transformations')
                for i = 1:M
                    hh = pi*(heading(i)-90)/180;
                    pp = pi*pitch(i)/180;
                    rr = pi*roll(i)/180;
                    
                    H = [cos(hh) sin(hh) 0; -sin(hh) cos(hh) 0; 0 0 1];
                    
                    % Make tilt matrix
                    P = [cos(pp) -sin(pp)*sin(rr) -cos(rr)*sin(pp);...
                        0             cos(rr)          -sin(rr);  ...
                        sin(pp) sin(rr)*cos(pp)  cos(pp)*cos(rr)];
                    % Make resulting transformation matrix depending on coord system
                    R = H*P*T;
                    %transform the data using rotated heading to correct for magnetic variation
                    % do each sample individually for now
                    V = [V1(i,:) V2(i,:) V3(i,:)]';
                    vel = R*V;
                    ADV.U(i,:) = vel(1);
                    ADV.V(i,:) = vel(2);
                    ADV.W(i,:) = vel(3);
                    if (i/100000)-floor(i/100000)==0
                        disp([num2str(i) ' samples of ' num2str(numel(heading)) ' transformed']),
                    end
                end
        end
end
ADV.Metadata.instmeta.coordsys = 'ENU'; %reassign

clearvars V1 V2 V3 vel V
disp(['Data rotated into Earth Coordinates using a magnetic variation of ',num2str(magvardeg)]);
end