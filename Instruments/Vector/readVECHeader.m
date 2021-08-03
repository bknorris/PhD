function Instmeta = readVECHeader(hdrfile)
%function to get instrument metadata about ADV deployment from header file

hdrFile = strcat(hdrfile,'.hdr');
hdr = fopen(hdrFile);
str='User setup';
while (~strncmp(str,'Hardware configuration',22));
     str=fgetl(hdr);
     if (strncmp(str,'Sampling rate',13))
         is=findstr(str,'Hz');
         Instmeta.VECSamplingRate = str2num(str(39:is-2));
     elseif (strncmp(str,'Nominal velocity range',22))
         is=findstr(str,'m/s');
         Instmeta.VECNominalVelocityRange = str2num(str(39:is-2));
     elseif (strncmp(str,'Burst interval',14))
         is=findstr(str,'sec');
         Instmeta.VECBurstInterval = str2num(str(39:is-2));
         if isempty(str(39:is-2))
             Instmeta.VECBurstInterval = 'CONTINUOUS';
         end
     elseif (strncmp(str,'Samples per burst',17))
         Instmeta.VECSamplesPerBurst = str2num(str(39:end));
         if strcmp(str(39:end),'N/A')
             Instmeta.VECSamplesPerBurst = 'N/A';
         end
     elseif (strncmp(str,'Sampling volume',15))
         is=findstr(str,'mm');
         Instmeta.VECSamplingVolume = str2num(str(39:is-2));
     elseif (strncmp(str,'Measurement Load',16))
         is=findstr(str,'%');
         Instmeta.VECMeasurementLoad = str2num(str(is-3:is-1));
     elseif (strncmp(str,'Transmit length',14))
         is=findstr(str,'mm');
         Instmeta.VECTransmitLength = str2num(str(39:is-2));
     elseif (strncmp(str,'Receive length',14))
         is=findstr(str,'m');
         Instmeta.VECReceiveLength = str2num(str(39:is-2));
     elseif (strncmp(str,'Analog output',13))
         Instmeta.VECAnalogOutputEnabled = (str(39:end));
     elseif (strncmp(str,'Analog input 1',14))
         Instmeta.VECAnalogInput1Sampling = (str(39:end));
     elseif (strncmp(str,'Analog input 2',14))
         Instmeta.VECAnalogInput2Sampling = (str(39:end));
     elseif (strncmp(str,'Power output',12))
         Instmeta.VECPowerOutputEnabled = (str(39:end));
     elseif (strncmp(str,'Output format',13))
         Instmeta.VECOutputFormat = (str(39:end));
     elseif (strncmp(str,'Velocity scaling',16))
         is=findstr(str,'mm');
         Instmeta.VECVelocityScaling = str2num(str(39:is-2));
     elseif (strncmp(str,'Powerlevel',10))
         Instmeta.VECPowerLevel = (str(39:end));
     elseif (strncmp(str,'Coordinate system',15))
         Instmeta.VECCoordinateSystem = str(39:end);
     elseif (strncmp(str,'Sound speed',11))
         Instmeta.VECSoundSpeed = str(39:end);
     elseif (strncmp(str,'Diagnostics measurements',17))
         Instmeta.VECDiagMeas = (str(39:end));
         %          elseif (strncmp(str,'Diagnostics - Interval',20))
         %              is=findstr(str,'sec');
         %              Instmeta.VECDiagInterval = str2num(str(39:is-2));
         %          elseif (strncmp(str,'Diagnostics - Number of samples',26))
         %              Instmeta.VECDiagNumberOfSamples = str2num(str(39:end));
         %          elseif (strncmp(str,'Diagnostics - Number of pings',26))
         %              Instmeta.VECDiagNumberOfPings = str2num(str(39:end));
     elseif (strncmp(str,'Salinity',8))
         is=findstr(str,'ppt');
         Instmeta.VECAssumedSalinity = str2num(str(39:is-2));
     elseif (strncmp(str,'Distance between pings',17))
         is=findstr(str,'m');
         Instmeta.VECistanceBetweenPings = str2num(str(39:is-2));
     elseif (strncmp(str,'Number of beams',15))
         Instmeta.VECNumberOfBeams = str2num(str(39));
     elseif (strncmp(str,'Software version',15))
         Instmeta.VECSoftwareVersion = num2str(str(39:42));
     elseif (strncmp(str,'Deployment name',15))
         Instmeta.VECeploymentName = str(39:end);
     elseif (strncmp(str,'Wrap mode',9))
         Instmeta.VECWrapMode = str(39:end);
     elseif (strncmp(str,'Deployment time',15))
         Instmeta.VECeploymentTime = str(39:end);
     elseif (strncmp(str,'Comments',8))
         Instmeta.VECComment1 = str(39:end);
         str=fgetl(hdr);
         if ~strncmp(str,'System1',7)
            Instmeta.VECComment2 = str(39:end);
         else
            str=fgetl(hdr);
            if ~strncmp(str,'System1',7)
                Instmeta.VECComment3 = str(39:end);
            else
                continue
            end
         end
     end
end
while (~strncmp(str,'Head configuration',18));
    str=fgetl(hdr);
    if (strncmp(str,'Serial number',13))
        Instmeta.VECInstrumentSerialNumber = str(39:end);
    elseif (strncmp(str,'Internal code version',21))
        Instmeta.VECInternalCodeVersion = num2str(str(39:end));
    elseif (strncmp(str,'Revision number',15))
        Instmeta.VECRevisionNumber = num2str(str(39:end));
    elseif (strncmp(str,'Recorder size',13))
        Instmeta.VECRevisionNumber = num2str(str(39:end));
    elseif (strncmp(str,'Firmware version',15))
        Instmeta.VECFirmwareVersion = num2str(str(39:end));
    end
end

while (~strncmp(str,'Data file format',16));
     str=fgetl(hdr);
     if strncmp(str,'Pressure sensor',15)
         Instmeta.PressureSensor = str(39:41);
     elseif strncmp(str,'Compass',7)
         Instmeta.Compass = str(39:41);
     elseif strncmp(str,'Tilt sensor',11)
         Instmeta.TiltSensor = str(39:41);
     elseif (strncmp(str,'Head frequency',12))
         is=findstr(str,'kHz');
         Instmeta.VECFrequency = str2num(str(39:is-2));
     elseif (strncmp(str,'Serial number',13))
         Instmeta.VECHeadSerialNumber = str(39:46);
     elseif (strncmp(str,'Transformation Matrix',15))
         Instmeta.VECTransformationMatrix = zeros(3,3);
         Instmeta.VECTransformationMatrix(1,:) = strread(str(39:end));
         str=fgetl(hdr);
         Instmeta.VECTransformationMatrix(2,:) = strread(str(39:end));
         str=fgetl(hdr);
         Instmeta.VECTransformationMatrix(3,:) = strread(str(39:end));
     elseif strncmp(str,'Pressure sensor calibration',27)
         Instmeta.PressureSensor = strread(str(39:end));
     end
end

fclose(hdr)
