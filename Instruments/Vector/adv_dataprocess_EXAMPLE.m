% This is an example metadata file for a Nortek Vector.
if ~exist('ADV','var')
    clear %a good idea
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

filename = 'vc101Sept27_3rdAttempt';
metadata.deployment_date = [2014 09 27 12 30 00]; %format [yyyy mm dd hh mm ss]
metadata.recovery_date = [2014 09 27 16 00 00]; %format [yyyy mm dd hh mm ss]
metadata.inst_lat = 9.491609; %latitude of the instrument
metadata.inst_lon = 106.244335; %longitude of the instrument
metadata.instname = 'VC101';
metadata.inst_type = 'Nortek Vector';            % type of instrument and instrument manufacturer
metadata.data_cmt = 'Mekong 2014 Control Experiment Vector VC101'; % any comment
metadata.vector_type = 'FIXED'; %or 'FLEXIBLE'
metadata.orientation = 'HORIZONTAL'; %'UP','DOWN', or 'HORIZONTAL' for the case
metadata.pressure_sensor_height =  450; %in mm
metadata.inst_hab = 450; %in mm
metadata.fixed_heading = 293; %compass measurement (in degrees) of the case
metadata.magdec = 0.26; %magnetic declination to be applied
% metadata.probe_height = []; %in mm

%If metadata.vector_type is 'FLEXIBLE', define properties of the PROBE
% metadata.case_rotation = 'UP'; %'DOWN' the direction of the Z axis
% metadata.probe_orientation = 'UP'; %'DOWN' - if different from the case rotation                    
% metadata.probe_heading = []; %compass measurement (in degrees) of the probe - if different from the case
% metadata.probe_pressure_dist = []; %distance between probe and case in mm (fixed is 214mm)
metadata.zero_pitch = 1; %'0' or '1'  override to zero pitch
metadata.zero_roll = 1; %'0' or '1'  override to zero roll
metadata.cutoff_amplitude = 90; %specify cutoff amplitude for correlations; sometimes the default (70) is too high

%If metadata.orientation is 'HORIZONTAL' swap beams to determine which beam
%is up. This is based on the deployment configuration, and needs to be
%determined before processing. This will be used in the conversion from XYZ
%to ENU coordinates only. Use a negative sign to flip the beam (*-1).
%Comment this section out if it will not be used.
metadata.b1 = 'b1'; %'format: 'b2' or 'b3', '-b2' or '-b3'
metadata.b2 = '-b2'; %'format: 'b1' or 'b3', '-b1' or '-b3'
metadata.b3 = 'b3'; %'format: 'b1' or 'b2', '-b1' or '-b2'

metdir = 'c:\Users\bkn5\Projects\Mekong_F2014\Data\Weather\';
diary(sprintf('run%s',datestr(now, 30)))

if 0
    ADV = readADVdata(filename,metadata);
end

%Downsample Pressure to 1Hz through interpolation to match temperature
%data, run pressure adjustment
if 0 
    Pres = interp1(ADV.Pres,1:ADV.Metadata.instmeta.samprate:length(ADV.Pres))';
    load([metdir 'metdata.mat'])
    metp = met.pressure.*0.0001; %convert Pa to dBar
    ind = find(met.datetime >= ADV.Sensor.Datetime(1) & met.datetime <= ADV.Sensor.Datetime(end));
    mett = met.datetime(ind);
    metp = metp(ind);
    ADV.Sensor.Pres = correctatmp(mett,metp,ADV.Sensor.Datetime,Pres,ADV.Sensor.Temp);
end

%Apply speed-of-sound correction
if 0
    t = ADV.Sensor.Temp;
    s = ADV.Metadata.instmeta.assumedsalinity;
    ss = ADV.Sensor.Soundspeed;
    [ADV.V1,ADV.V2,ADV.V3] = adjsound(t,s,ss,ADV.V1,ADV.V2,ADV.V3);
end

%Rotate velocities to ENU, optionally run despike routine
if 0
   ADV = rotateADV(ADV,metadata);
   ADV.Metadata.instmeta.coordsys = 'ENU'; %reassign
end

%run despike program
if 1
    [U,V,W] = despikeADV(ADV,metadata);
end

%assign variables
if 1 
    ADV.U = U; ADV.V = V; ADV.W = W;
    clearvars U V W
end

%Check Amplitudes, calculate Statistics
if 1
   ADV = ADVampcheck(ADV,metadata); 
end

%Save data file
if 0
    dd = datestr(ADV.Metadata.deployment_date,'ddmmyy');
    fname = [ADV.Metadata.instname '_' dd];
    save(fname,'ADV') %save final file; f for "final"
    disp(['file saved as ' fname '.mat'])
    clearvars -except ADV metadata filename
end

if 1
    plotDespikecmp
end

diary off