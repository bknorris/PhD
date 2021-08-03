function [stereoOut f_coeff] = dopplerDF(input, v, y0, Fs)
%dopplerDF - a process of file with Doppler Effect
%This function applies the Doppler Effect to the input vector due to
%the moving sound source.
%
% Syntax:  output f_coeff = doppler_effect(input, speed, mindist, Fs)
%
% Inputs:
%    input    - input vector to be processed by the function
%    speed    - Doopler Effect depends on speed of the moving sound source,
%               that is why the speed needs to be specified (metres per second)
%    mindist  - the distance represented by perpendicular vector of the
%               observer to the line of motion of the sound source (metres)
%    Fs       - sampling frequency of the sound vector 'input'
%
% Outputs:
%    output   - a vector storing the processed file
%    f_coeff  - a vector of the frequency change coefficients
%

%% generating vector of frequency changes due to doppler effect
input=input';
N = length(input); % number of samples
t = N/Fs; % duration of input
x0 = v*t; % distance covered by source
x = linspace(-x0/2, x0/2, N); % vector of horizontal displacement
alpha = atan2(y0,-x); % angle between source motion and line from source to observer

vdoppler = v.*cos(alpha); % radial velocity changes due to the angle between source and observer
c = 343; % sound velocity
v0 = 0; % the velocity of the observer is = 0 metres per second

f_coeff = (c + vdoppler)./(c - v0); % vector of frequency change coefficients

%% Applying doppler effect to resample wave file

chunksize = 500; % length of one chunk
overlap = 0; % percentage overlap
numChunks = ceil(N/(chunksize*(1-(overlap/100)))); % number of chunks
output = 0; % initialise output
for i = 0:numChunks-1
     istart = (i * chunksize + 1)-(chunksize*i*(overlap/100)); 
     iend = istart + chunksize;
     if iend > N; 
         iend = N;
     end
     fcoeff_av = floor(mean(f_coeff(istart:iend))*100)/100;
     output = [output resample(input(istart:iend), 100, floor(fcoeff_av*100))];
end
output=output./max(abs(output)); % Normalise

%% 2 channel output 

iad = 0.17; % Interaural distance in metres
iad = 4*iad; % note: increased by factor of 4 to make effect more obvious

% Left Ear
distanceL = abs((x+iad/2)./(cos(alpha))); % vector of distance to sound source
attenL = 1./(1+distanceL); % vector of attenuation factors based on distance
% note: 1 is added to the denomenator to avoid divide by zero errors

% Right Ear
distanceR = abs((x-iad/2)./(cos(alpha))); % vector of distance to sound source
attenR = 1./(1+distanceR); % vector of attenuation factors based on distance

% use fixes so that matrix dimensions match for multiplication
fix1 = floor((length(output)-length(attenL))/2);
fix2 = ceil((length(output)-length(attenL))/2);
Fix1 = zeros(1, fix1);
Fix2 = zeros(1, fix2);

% Apply attenduation  factors to output (doppler frequency shifted)
outputL=output.*[Fix1, attenL, Fix2];
outputR=output.*[Fix1, attenR, Fix2];

% Combine left and right for stereo output
stereoOut = [outputL; outputR]';

% Normalise output
stereoOut = stereoOut./max(max(abs(stereoOut)));

end


