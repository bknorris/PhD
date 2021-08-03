close all; clear;

% generate time array of XXX days
dt = 10*60; % = 10 min sampling
N = 4320; % 10 min sampling for 30 days yields N = 4320 samples 
N/2;
t = 0:dt:N*dt; % sec

% M2 tide
T_m2 = 12.42*60*60;
A_m2 = 1;
function_m2 = A_m2*cos(2*pi*t/T_m2);

% M4 tide
% T_m4 = 6.21*60*60;
% A_m4 = 1;
% function_m4 = 0*A_m4*cos(2*pi*t/T_m4 - 0*2*pi/360*314.65);

% S2 tide
T_s2 = 12*60*60;
A_s2 = 1;
function_s2 = A_s2*cos(2*pi*t/T_s2);

% MSF tide
% T_msf = 354.3712*60*60;
% A_msf = 0.5;
% function_msf = A_msf*cos(2*pi*t/T_msf- 0*2*pi/360*310.8);

combined = function_m2 + function_s2;
plot(t,combined)
grid on

% generate frequencies
f0 = 1/(N*dt);
fN = 1/(2*dt);
df = f0;
f = f0:df:fN;
fN/f0;

% get energy
Xn=fft(combined,N);
Pxx = Xn.*conj(Xn)/N;
% get rid of negative freq parts of Pxx
Pxx(N/2+1:end) = [];
% and double power of pos freq parts
Pxx(2:end) = 2*Pxx(2:end);

figure


hold on
for i = 1:length(f)
plot([f(i) f(i)],[0 max(Pxx)],'k')
end

plot(f,Pxx,'m','linewidth',[3])
grid
xlabel('frequency (Hz)')
ylabel('Energy (m^2/Hz  or m^2 s)')

