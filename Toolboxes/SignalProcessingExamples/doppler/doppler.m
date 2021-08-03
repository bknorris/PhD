%   NymFys - Example 2.1 - Doppler effect
%   31 December 2012

% Call the FFT wrapper function
[P1,P2,f]=fftwrapper(true);

% Plot power spectrum
plot(f,P1./max(P1),f,P2./max(P2));
xlabel('f (Hz)');
ylabel('P/P_{max}')
legend('Sample 1','Sample 2','Location','NorthWest');

% Find frequency peaks (this works for this example, but NOT in general)
f1=f(P1==max(P1))
f2=f(P2==max(P2))

c=340.29; % Speed of sound [m/s]
v=(f1-f2)/(f1+f2)*c; % Speed of train [m/s]
v=v*3.6; % Speed of of train [km/h]

fprintf('f1 : %3.1f Hz\n',f1);
fprintf('f2 : %3.1f Hz\n',f2);
fprintf('v  : %3.2f km/h\n',v);