function [Instmeta,SENFIELDS] = readNortekGenericHeader(hdrfile)
%function to get instrument metadata about AWAC deployment from header file

hdrFile = strcat(hdrfile,'.hdr');
hdr = fopen(hdrFile);
str='User setup';
while (~strncmp(str,'Hardware configuration',22));
     str=fgetl(hdr);
     if (strfind(str,'Profile interval'))
        is=findstr(str,'sec');
        Instmeta.ProfileInterval = str2num(str(39:is-2));
     elseif (strncmp(str,'Sampling rate',13))
        is=findstr(str,'Hz');
        Instmeta.SamplingRate = str2num(str(39:is-2));
     elseif (strncmp(str,'Burst interval',14))
        is=findstr(str,'sec');
     elseif (strncmp(str,'Measurement/Burst interval',26))
        is=findstr(str,'sec');
        Instmeta.BurstInterval = str2num(str(39:is-2));
     elseif (strncmp(str,'Samples per burst',17))
        Instmeta.SamplesPerBurst = str2num(str(39:end));
     end
end

%grab the field numbers for Pitch, Roll and Heading, as it can vary from
%one instrument to another
while ~feof(hdr)
    str = fgetl(hdr);
    if any(strfind(str,'.sen'))
        while (~strcmp(str,'---------------------------------------------------------------------'));
            str = fgetl(hdr);
            if (strfind(str,'Heading'))
                X = textscan(str,'%4c%33c\n');
                SENFIELDS.Heading = str2num(X{1});
            elseif (strfind(str,'Pitch'))
                X = textscan(str,'%4c%33c\n');
                SENFIELDS.Pitch = str2num(X{1});
            elseif (strfind(str,'Roll'))
                X = textscan(str,'%4c%33c\n');
                SENFIELDS.Roll = str2num(X{1});
            else
            end
        end
    end
end

fclose(hdr);