%now do rotations from either Beam or XYZ into ENU
disp('Performing transformations to Earth Coordinates...')
%transformation matrix for S/N 
T = ncread(rawfile,'TransMatrix');
% it seems that some old aquadopp data has a different transmatrix
% old stuff is integer
if isinteger(T(1,1)), % integer that needs scaling
    fprintf('The TransMatrix has been scaled by T=T/4096\n')
    T = T/4096;   % Scale the transformation matrix correctly to floating point numbers
end

%
VEL.U = zeros(M,N);VEL.V = zeros(M,N);VEL.W = zeros(M,N);

%correct heading for magnetic variation
% use whichever of cdf.magnetic_variation(:) or cdf.magnetic_variation_at_site(:)
% problem with this is that it doesn't throw an error if it doesn't exist,
% just an empty matrix
try
    magvardeg = Gatts.magnetic_variation_at_site(:);
catch iter
end
if isempty(magvardeg)
    magvardeg = Gatts.magnetic_variation(:);
end
if isempty(magvardeg)
    magvardeg = 0;
end

% see if the user wants to use a fixed heading value
if isfield(metadata,'fixed_heading')
    heading = ones(N,1)*metadata.fixed_heading;
end

fprintf('Heading is being rotated by %f\n',magvardeg);
heading = heading + magvardeg;
heading( heading >= 360 ) = heading(heading >= 360) - 360;
heading( heading < 0 ) = heading(heading < 0) + 360;


magvar = magvardeg(1)*pi/180;

switch Gatts.VECCoordinateSystem
    case 'ENU'
        disp('Data are already in Earth Coordinates, applying magnetic variation')
        VEL.U =  vel1*cos(magvar) + vel2*sin(magvar);
        VEL.V = -vel1*sin(magvar) + vel2*cos(magvar);
        VEL.W = vel3;
    case 'XYZ'
        disp('Data are in XYZ coordinates - applying rotations')
        for i = 1:N
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
            for j = 1:spb
                V = [vel1(j,i) vel2(j,i) vel3(j,i)]';
                vel = R*V;
                VEL.U(j,i) = vel(1);
                VEL.V(j,i) = vel(2);
                VEL.W(j,i) = vel(3);
            end
        end
            %below is the old way of doing it
%             VEL.U =  U*cos(magvar) + V*sin(magvar);
%             VEL.V = -U*sin(magvar) + V*cos(magvar);
    case 'BEAM'
        disp('Data are in BEAM coordinates - applying beam transformations')
        for i = 1:N
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
            for j = 1:spb
                V = [vel1(j,i) vel2(j,i) vel3(j,i)]';
                vel = R*V;
                VEL.U(j,i) = vel(1);
                VEL.V(j,i) = vel(2);
                VEL.W(j,i) = vel(3);
            end
        end
            %below is the old way of doing it
%             VEL.U =  U*cos(magvar) + V*sin(magvar);
%             VEL.V = -U*sin(magvar) + V*cos(magvar);
end
rotcomment = ['Data rotated into Earth Coordinates using a magnetic variation of ',num2str(magvardeg),' by vec2ncBurstNative. '];
disp(rotcomment)
history = [rotcomment ' | ' history];


%check amplitudes - Nortek lists the noise floor as 70 counts for the
%vector
VEL.AGC1 = AMP1;VEL.AGC2 = AMP2;VEL.AGC3 = AMP3;
VEL.pitch = pitch;VEL.roll = roll;VEL.heading = heading;VEL.pressure = press;VEL.temperature = temp;

disp('Screening data by predetermined Amplitude cutoff for this frequency system....')
if isfield(metadata,'amplitude_cutoff')
    cutoff = metadata.amplitude_cutoff;
    disp(['User selected an amplitude cutoff of ',num2str(cutoff),' for screening data'])
else
    cutoff = 70;
end
ampcomment = ['Amplitude of ',num2str(cutoff),' counts used for data screening'];
history = [ampcomment ' | ' history];

ind = find(AMP1<cutoff | AMP2<cutoff | AMP3<cutoff);
VEL.U(ind) = NaN;
VEL.V(ind) = NaN;
VEL.W(ind) = NaN;

% calculate the statistics; for time, we will use the middle of the burst
flds = fieldnames(VEL);
for i = 1:length(flds)
    fld = flds{i};
    if length(VEL.(fld))>1; %i.e. only work on timeseries
        STATS.(fld) = nanmean(VEL.(fld));
        if any(strcmp(fld,{'U';'V';'W';'pressure'}))
            stdnm = [fld '_STD'];STATS.(stdnm) = nanstd(VEL.(fld))';
            maxnm = [fld '_MAX'];STATS.(maxnm) = nanmax(VEL.(fld))';
            minnm = [fld '_MIN'];STATS.(minnm) = nanmin(VEL.(fld))';
        end
    end
end
STATS.time = floor(STATS.jd);STATS.time2 = (STATS.jd-floor(STATS.jd))*86400000;

clear vel1 vel2 vel3 pitch roll heading temp press AMP1 AMP2 AMP3

% % put in fill values
% N = length(VEL.jd);flds = fieldnames(VEL);
% for i = 1:length(VEL)
%     fld = flds{i};
%     if length(VEL.(fld) == N); %i.e. only work on timeseries
%         if strcmp(fld,{'time';'time2'}),
%             continue
%         else
%             nanind = isnan(VEL.(fld));VEL.(fld)(nanind) = doublefill;
%         end
%     end
% end