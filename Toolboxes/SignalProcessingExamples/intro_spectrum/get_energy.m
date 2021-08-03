function [S,amp_squared,f]=get_energy(A,B,N,dt)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

%% generate frequencies
f0 = 1/(N*dt);  % fundamental frequnecy
fN = 1/(2*dt);  % Nyquist frequency
df = f0;        % frequenct increment
f = f0:df:fN;   % frequency array

%% get "energy"
%figure
%subplot(2,1,1)
amp_squared = A.^2 + B.^2;
%bar(f,E,'k')
%grid
%xlabel('f (Hz)')
%ylabel('Amplitude Squared (units^2)')


%% get energy density
S = 1/2*amp_squared/df;
%subplot(2,1,2)
%plot(f,S,'k');
%grid
%xlabel('f (Hz)')
%ylabel('Energy Density (units^2/Hz or units^2*s)')
%disp(['Total signal energy is ',num2str(trapz(f,S)), ' units^2'])




end

 