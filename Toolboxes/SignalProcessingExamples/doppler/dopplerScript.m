clear all
clc
%% Read in wave file and input parameters
sample = 1; % choice between sample 1, 2 or 3 using the swith case
switch sample
    case 1
        [inwave fs nbits]=wavread('Air_Horn.wav');
        speed = 20; % speed of sound source in metres per second
        mindist = 2; % minimum distance of sound source from observer in metres
    case 2
        [inwave fs nbits]=wavread('Siren.wav');
        speed = 30; % speed of sound source in metres per second
        mindist = 5; % minimum distance of sound source from observer in metres
    case 3
        [inwave fs nbits]=wavread('Helicopter.wav');
        speed = 50; % speed of sound source in metres per second
        mindist = 50; % minimum distance of sound source from observer in metres
end

%% Function call
[outwave f] = dopplerDF(inwave, speed, mindist, fs); 

%% Play sound
sound(outwave, fs);

%% Save output as wave file
wavwrite(outwave, fs, 'dopplerOutput.wav'); 