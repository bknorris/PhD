function Instmeta = readVEChdr(filename)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%A function to get instrument metadata about ADV deployment from header file,
%and user set up from the sen file. This function creates the structure adv
%and populates it with values from both the header and sen files, then
%loads the velocity data and saves it to this structure.
%
%  Contains adaptations of readNortekGenericHeader.m and readVECHeader.m
%  Written by Kurt J Rosenberger
%  krosenberger@usgs.gov
%  USGS Pacific Coastal Marine Science Center
%
%  This script was compiled by Benjamin K Norris, 2015
%  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('Reading the header file')
hdrFile = strcat(filename,'.hdr');
hdr = fopen(hdrFile);
top = textscan(hdr,'%s',9,'delimiter','\n'); %read the first 10 rows of the hdr file
%check for checksum errors
str = 'Number of velocity checksum errors';
checksum = strncmpi(top{1},str,34);
num = regexprep(top{1}(checksum),str,'');
if str2num(num{1}) ~= 0
    disp('There are checksum errors')
else
    disp('There are no checksum errors')
end
clearvars top num
str='User setup';
while (~strncmp(str,'Hardware configuration',22));
     str=fgetl(hdr);
     if (strncmp(str,'Sampling rate',13))
         is=findstr(str,'Hz');
         Instmeta.samprate = str2num(str(39:is-2));
     elseif (strncmp(str,'Nominal velocity range',22))
         is=findstr(str,'m/s');
         Instmeta.nominalvelrange = str2num(str(39:is-2));
     elseif (strncmp(str,'Burst interval',14))
         is=findstr(str,'sec');
         Instmeta.burstint = str2num(str(39:is-2));
         if isempty(str(39:is-2))
             Instmeta.burstint = 'CONTINUOUS';
             disp('Instrument recorded in CONTINUOUS mode')
         else
             disp('Instrument recorded in BURST mode')
         end
     elseif (strncmp(str,'Samples per burst',17))
         Instmeta.spb = str2num(str(39:end));
         if strcmp(str(39:end),'N/A')
             Instmeta.spb = 'N/A';
         end
     elseif (strncmp(str,'Sampling volume',15))
         is=findstr(str,'mm');
         Instmeta.sampvol = str2num(str(39:is-2));
     elseif (strncmp(str,'Measurement Load',16))
         is=findstr(str,'%');
         Instmeta.measload = str2num(str(is-3:is-1));
     elseif (strncmp(str,'Transmit length',14))
         is=findstr(str,'mm');
         Instmeta.trasmlength = str2num(str(39:is-2));
     elseif (strncmp(str,'Receive length',14))
         is=findstr(str,'m');
         Instmeta.recievelength = str2num(str(39:is-2));
     elseif (strncmp(str,'Analog output',13))
         Instmeta.analogoutputenabled = (str(39:end));
     elseif (strncmp(str,'Analog input 1',14))
         Instmeta.analoginput1sampling = (str(39:end));
     elseif (strncmp(str,'Analog input 2',14))
         Instmeta.analoginput2sampling = (str(39:end));
     elseif (strncmp(str,'Power output',12))
         Instmeta.poweroutputenabled = (str(39:end));
     elseif (strncmp(str,'Output format',13))
         Instmeta.outputformat = (str(39:end));
     elseif (strncmp(str,'Velocity scaling',16))
         is=findstr(str,'mm');
         Instmeta.velocityscaling = str2num(str(39:is-2));
     elseif (strncmp(str,'Powerlevel',10))
         Instmeta.powerlevel = (str(39:end));
     elseif (strncmp(str,'Coordinate system',15))
         Instmeta.coordsys = str(39:end);
     elseif (strncmp(str,'Sound speed',11))
         Instmeta.soundspd = str(39:end);
     elseif (strncmp(str,'Diagnostics measurements',17))
         Instmeta.diagnosticmsg = (str(39:end));
         %          elseif (strncmp(str,'Diagnostics - Interval',20))
         %              is=findstr(str,'sec');
         %              Instmeta.VECDiagInterval = str2num(str(39:is-2));
         %          elseif (strncmp(str,'Diagnostics - Number of samples',26))
         %              Instmeta.VECDiagNumberOfSamples = str2num(str(39:end));
         %          elseif (strncmp(str,'Diagnostics - Number of pings',26))
         %              Instmeta.VECDiagNumberOfPings = str2num(str(39:end));
     elseif (strncmp(str,'Salinity',8))
         is=findstr(str,'ppt');
         Instmeta.assumedsalinity = str2num(str(39:is-2));
     elseif (strncmp(str,'Distance between pings',17))
         is=findstr(str,'m');
         Instmeta.distancebetweenpings = str2num(str(39:is-2));
     elseif (strncmp(str,'Number of beams',15))
         Instmeta.numofbeams = str2num(str(39));
     elseif (strncmp(str,'Software version',15))
         Instmeta.softwareversion = num2str(str(39:42));
     elseif (strncmp(str,'Deployment name',15))
         Instmeta.depname = str(39:end);
     elseif (strncmp(str,'Wrap mode',9))
         Instmeta.wrapmode = str(39:end);
     elseif (strncmp(str,'Deployment time',15))
         Instmeta.deptime = str(39:end);
     elseif (strncmp(str,'Comments',8))
         Instmeta.comment1 = str(39:end);
         str=fgetl(hdr);
         if ~strncmp(str,'System1',7)
            Instmeta.comment2 = str(39:end);
         else
            str=fgetl(hdr);
            if ~strncmp(str,'System1',7)
                Instmeta.comment3 = str(39:end);
            else
                continue
            end
         end
     end
end
while (~strncmp(str,'Head configuration',18));
    str=fgetl(hdr);
    if (strncmp(str,'Serial number',13))
        Instmeta.serialno = str(39:end);
    elseif (strncmp(str,'Internal code version',21))
        Instmeta.internalcodeversion = num2str(str(39:end));
    elseif (strncmp(str,'Revision number',15))
        Instmeta.revisionno = num2str(str(39:end));
    elseif (strncmp(str,'Recorder size',13))
        Instmeta.recsize = num2str(str(39:end));
    elseif (strncmp(str,'Firmware version',15))
        Instmeta.firmwareversion = num2str(str(39:end));
    end
end

while (~strncmp(str,'Data file format',16));
     str=fgetl(hdr);
     if strncmp(str,'Pressure sensor',15)
         Instmeta.pressen = str(39:41);
     elseif strncmp(str,'Compass',7)
         Instmeta.compass = str(39:41);
     elseif strncmp(str,'Tilt sensor',11)
         Instmeta.tiltsen = str(39:41);
     elseif (strncmp(str,'Head frequency',12))
         is=findstr(str,'kHz');
         Instmeta.frequency = str2num(str(39:is-2));
     elseif (strncmp(str,'Serial number',13))
         Instmeta.headserialno = str(39:46);
     elseif (strncmp(str,'Transformation Matrix',15))
         Instmeta.transfm = zeros(3,3);
         Instmeta.transfm(1,:) = strread(str(39:end));
         str=fgetl(hdr);
         Instmeta.transfm(2,:) = strread(str(39:end));
         str=fgetl(hdr);
         Instmeta.transfm(3,:) = strread(str(39:end));
     elseif strncmp(str,'Pressure sensor calibration',27)
         Instmeta.pressencal = strread(str(39:end));
     end
end

%grab the field numbers for the sensor data, as they can vary from
%one instrument to another
while ~feof(hdr)
    str = fgetl(hdr);
    if any(strfind(str,'.sen'))
        while (~strcmp(str,'---------------------------------------------------------------------'));
            str = fgetl(hdr);
            if (strfind(str,'Day'))
                X = textscan(str,'%4c%33c\n');
                SENFIELDS.Day = str2num(X{1});
            elseif (strfind(str,'Month'))
                X = textscan(str,'%4c%33c\n');
                SENFIELDS.Month = str2num(X{1});
            elseif (strfind(str,'Year'))
                X = textscan(str,'%4c%33c\n');
                SENFIELDS.Year = str2num(X{1});
            elseif (strfind(str,'Hour'))
                X = textscan(str,'%4c%33c\n');
                SENFIELDS.Hour = str2num(X{1});
            elseif (strfind(str,'Minute'))
                X = textscan(str,'%4c%33c\n');
                SENFIELDS.Minute = str2num(X{1});
            elseif (strfind(str,'Second'))
                X = textscan(str,'%4c%33c\n');
                SENFIELDS.Second = str2num(X{1});
            elseif (strfind(str,'Error code'))
                X = textscan(str,'%4c%33c\n');
                SENFIELDS.Errorcode = str2num(X{1});
            elseif (strfind(str,'Status code'))
                X = textscan(str,'%4c%33c\n');
                SENFIELDS.Statuscode = str2num(X{1});
            elseif (strfind(str,'Battery voltage'))
                X = textscan(str,'%4c%33c\n');
                SENFIELDS.Battvolt = str2num(X{1});
            elseif (strfind(str,'Soundspeed'))
                X = textscan(str,'%4c%33c\n');
                SENFIELDS.Soundspeed = str2num(X{1});
            elseif (strfind(str,'Heading'))
                X = textscan(str,'%4c%33c\n');
                SENFIELDS.Heading = str2num(X{1});
            elseif (strfind(str,'Pitch'))
                X = textscan(str,'%4c%33c\n');
                SENFIELDS.Pitch = str2num(X{1});
            elseif (strfind(str,'Roll'))
                X = textscan(str,'%4c%33c\n');
                SENFIELDS.Roll = str2num(X{1});
            elseif (strfind(str,'Temperature'))
                X = textscan(str,'%4c%33c\n');
                SENFIELDS.Temp = str2num(X{1});
            elseif (strfind(str,'Analog input'))
                X = textscan(str,'%4c%33c\n');
                SENFIELDS.Analoginput = str2num(X{1});
            elseif (strfind(str,'Checksum'))
                X = textscan(str,'%4c%33c\n');
                SENFIELDS.Checksum = str2num(X{1});
            else
            end
        end
    end
end
Instmeta.senflds = SENFIELDS;
fclose(hdr);

end