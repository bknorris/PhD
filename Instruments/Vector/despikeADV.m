function [V1,V2,V3] = despikeADV(ADV,metadata,p)
%This process is completed firstly by NaN-ing out
%velocities corresponding to a threshold of low correlations (<40%), then
%runs a despike routine based on the method outlined by Goring & Nikora
%(2002), modified by Nobuhito Mori (2005). For reference: 
%http://www.mathworks.com/matlabcentral/fileexchange/15361-despiking
%
% Inputs: ADV - the structure ADV containing metadata fields and velocity
%         data
%         metadata - additional metadata fields
%         plot     - optional, set to 1 if user wants to see the despike
%                    comparison (otherwise set to []). 
%
% Outputs: V1,V2,V3: despiked velocities
%
% Dependencies: this function calls: 
%               func_despike_phasespace3d_3var.m
%               func_despike_phasespace.m.
%
%   This script was compiled by Benjamin K Norris, 2015
%   University of Waikato, New Zealand
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('Running signal processing - despiking routine')
disp('DespikeADV needs to rotate XYZ or ENU coordinates into BEAM')

%set up the transformation matrix, convert xyz to beam coordinates
T = ADV.Metadata.instmeta.transfm;
samprate = ADV.Metadata.instmeta.samprate;
n = ADV.Metadata.nsamp*samprate;
beam = zeros(n,3);

%define h,p,r same as rotateADV; rotate to beam to perform despiking
if isfield(metadata,'fixed_heading') %user defined fixed heading value
    disp(['User requests a fixed heading of ' num2str(metadata.fixed_heading) ' degrees']);
    heading = ones(n,1)*metadata.fixed_heading;
else
    heading = interp(ADV.Sensor.Heading,samprate);
end
if metadata.zero_pitch %if user wants zero'd pitch/roll
    pitch = zeros(n,1);
else
    pitch = interp(ADV.Sensor.Pitch,samprate);
end
if metadata.zero_roll
    roll = zeros(n,1);
else
    roll = interp(ADV.Sensor.Roll,samprate);
end

switch ADV.Metadata.instmeta.coordsys
    case 'ENU'
        disp('Data are in Earth Coordinates - rotating to BEAM')
        for i = 1:n
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
            V = [ADV.U(i,:) ADV.V(i,:) ADV.W(i,:)]';
            xyz = (T*inv(R)*V); %convert from ENU to XYZ
            beam(i,:) = (inv(T)*xyz)'; %#ok<*MINV> convert from XYZ to BEAM
        end
        B1 = beam(:,1);B2 = beam(:,2);B3 = beam(:,3);
        disp('Velocities have been converted from ENU to BEAM')
    case 'XYZ'
        disp('Data are in XYZ coordinates - rotating to BEAM')
        for i = 1:n
            V = [ADV.V1(i,:) ADV.V2(i,:) ADV.V3(i,:)]';
            beam(i,:) = (inv(T)*V)'; %#ok<*MINV>
        end
        B1 = beam(:,1);B2 = beam(:,2);B3 = beam(:,3);
        disp('Velocities have been converted from XYZ to BEAM')
    case 'BEAM'
        disp('Data are all ready in BEAM coordinates')
end
clearvars beam

[m,n] = size(ADV.Cor1);

%filter out low correlations
cutoff = metadata.cutoff_amplitude;
ind = find(ADV.Cor1<cutoff | ADV.Cor2<cutoff | ADV.Cor3<cutoff);
B1(ind) = NaN;
B2(ind) = NaN;
B3(ind) = NaN;
%filter low SNR values
acutoff = 15;
ind = find(ADV.Snr1<acutoff | ADV.Snr2<acutoff | ADV.Snr3<acutoff);
B1(ind) = NaN;
B2(ind) = NaN;
B3(ind) = NaN;
beam1 = B1;beam2 = B2;beam3 = B3;
clear B1 B2 B3  

%report number of bad values found by script
invalid = find(isnan(beam1) | isnan(beam2) | isnan(beam3));
disp(['Eliminated ' num2str(length(invalid)) ' bad values in velocity data']);
disp(['Velocity measurements corresponding to correlations <' num2str(metadata.cutoff_amplitude) '% have been removed'])
disp('Velocity measurements corresponding to SNR values <15dB have been removed')

%bridge gaps in timeseries
tic
B1 = zeros(m,n);B2 = zeros(m,n);B3 = zeros(m,n);
disp('Bridging gaps in velocity time series')
samprate = ADV.Metadata.instmeta.samprate; %in Hz (samples/second)
nlin = samprate*30; %# of samples in 30sec
nmaxbr = nlin*8; %# of samples in 4 min
burst = 10*samprate*60; %# of samples in 10 minutes
gap1 = cmgidgaps(beam1);
gap2 = cmgidgaps(beam2);
gap3 = cmgidgaps(beam3);
if any([gap1 gap2 gap3] > 1)
    %crop beams to the same size;
    beam1 = beam1(1:m,1:n);
    beam2 = beam2(1:m,1:n);
    beam3 = beam3(1:m,1:n);
    id = [1:burst:m m];
    for i = 1:length(id)-1
        B1(id(i):id(i+1),:) = cmgbridge(beam1(id(i):id(i+1),:),nlin,nmaxbr,burst);
        B2(id(i):id(i+1),:) = cmgbridge(beam2(id(i):id(i+1),:),nlin,nmaxbr,burst);
        B3(id(i):id(i+1),:) = cmgbridge(beam3(id(i):id(i+1),:),nlin,nmaxbr,burst);
    end
end
%sometimes cmgbridge creates really large values (e.g. >100), find
%these and set them equal to zero
B1(B1>100) = 0;B1(B1<-100) = 0;
B2(B2>100) = 0;B2(B2<-100) = 0;
B3(B3>100) = 0;B3(B3<-100) = 0;

%Sometimes NaNs are missed... fill with zeros
invalid = find(isnan(B1) | isnan(B2) | isnan(B3));
B1(invalid) = 0;B2(invalid) = 0;B3(invalid) = 0;

clear beam1 beam2 beam3 
disp(['Gaps in data bridged in ' num2str(toc/60) ' minutes'])
%run despike routine
tic
disp('Running Despiking Routine')
[V1,~] = func_despike_phasespace3d(B1,[],2);
[V2,~] = func_despike_phasespace3d(B2,[],2);
[V3,~] = func_despike_phasespace3d(B3,[],2);
disp(['Spike noise interpolated in ' num2str(toc/60) ' minutes'])

if p == 1
    figure
    ax1 = subplot(311);
    plot(B1,'k')
    hold on
    plot(V1,'r')
    legend('No Despike','With Despike')
    title('V1')
    xlabel('N samples')
    ylabel('m/s')
    ax2 = subplot(312);
    plot(B2,'k')
    hold on
    plot(V2,'r')
    legend('No Despike','With Despike')
    title('V2')
    xlabel('N samples')
    ylabel('m/s')
    ax3 = subplot(313);
    plot(B3,'k')
    hold on
    plot(V3,'r')
    legend('No Despike','With Despike')
    title('V3')
    xlabel('N samples')
    ylabel('m/s')
    linkaxes([ax1 ax2 ax3],'xy')
    set([ax1 ax2 ax3],'Ylim',[-2 2])
end

%rotate BEAM velocities back to XYZ
disp('Re-rotating BEAM velocities back to XYZ...')
xyz = zeros(length(V1),3);
for i = 1:length(V1)
    V = [V1(i,:) V2(i,:) V3(i,:)]';
    xyz(i,:) = T*V;
end
V1 = xyz(:,1);V2 = xyz(:,2);V3 = xyz(:,3);
end