%generate the noisy signal which will be filtered
x= cos(2*pi*12*[0:0.001:1.23]);
x(end) = [];
[b a] = butter(2,[0.6 0.7],'bandpass');
filtered_noise = filter(b,a,randn(1, length(x)*2));
x = (x + 0.5*filtered_noise(500:500+length(x)-1))/length(x)*2;
% This script is available at http://dadorran.wordpress.com search for
% filtering matlab demo
plot(x)
title('Noisy signal')
xlabel('Samples');
ylabel('Amplitude')
 
%plot magnitude spectrum of the signal
X_mags = abs(fft(x));
figure(10)
plot(X_mags)
xlabel('DFT Bins')
ylabel('Magnitude')
 
%plot first half of DFT (normalised frequency)
num_bins = length(X_mags);
plot([0:1/(num_bins/2 -1):1], X_mags(1:num_bins/2))
xlabel('Normalised frequency (\pi rads/sample)')
ylabel('Magnitude')
 
%design a second order filter using a butterworth design technique 
[b a] = butter(2, 0.3, 'low')
 
%Show the signal flow diagram of the filter/system
%this function available from http://dadorran.wordpress.com
%create_signal_flow(b,a) 
 
%plot the frequency response (normalised frequency)
H = freqz(b,a, floor(num_bins/2));
hold on
plot([0:1/(num_bins/2 -1):1], abs(H),'r');
 
%filter the signal using the b and a coefficients obtained from
%the butter filter design function
x_filtered = filter(b,a,x);
 
%plot the filtered signal
figure(2)
plot(x_filtered,'r')
title('Filtered Signal - Using Second Order Butterworth')
xlabel('Samples');
ylabel('Amplitude')
 
%Redesign the filter using a higher order filter
[b2 a2] = butter(20, 0.3, 'low')
 
%Plot the magnitude spectrum and compare with lower order filter
H2 = freqz(b2,a2, floor(num_bins/2));
figure(10)
hold on
plot([0:1/(num_bins/2 -1):1], abs(H2),'g');
 
%filter the noisy signal and plot the result
x_filtered2 = filter(b2,a2,x);
figure(3)
plot(x_filtered2,'g')
title('Filtered Signal - Using 20th Order Butterworth')
xlabel('Samples');
ylabel('Amplitude')
 
%Use a band reject filter in place of a low pass filter
[b_stop a_stop] = butter(20, [0.5 0.8], 'stop');
 
%plot the magnitude spectrum
H_stopband = freqz(b_stop,a_stop, floor(num_bins/2));
figure(10)
hold on
plot([0:1/(num_bins/2 -1):1], abs(H_stopband),'c');
 
%plot filtered signal
x_filtered_stop = filter(b_stop,a_stop,x);
figure(4);
plot(x_filtered_stop,'c')
title('Filtered Signal - Using Stopband')
xlabel('Samples');
ylabel('Amplitude')
 
%Use matlabs built-in buttord function to get the optimum order to meet a specification
[N Wn] = buttord(0.1, 0.5, 5, 40)
 
%use the N and Wn values obtained above to design the filter in the usual way
[b3 a3] = butter(N, Wn, 'low');
 
%plot the magnitude spectrum
H3 = freqz(b3,a3, floor(num_bins/2));
figure(10);
hold on
plot([0:1/(num_bins/2 -1):1], abs(H2),'k');
figure(10)
 
%filter the signal and plot the ouput of the filter
x_filtered3 = filter(b3,a3,x);
figure(5);
plot(x_filtered3,'k')
title(['Filtered Signal - Using ' num2str(N) ' th Order Butterworth'])
xlabel('Samples');
ylabel('Amplitude')
% script available from http://dadorran.wordpress.com
% comparison with other filter design techniques (chebyshev and elliptical)
[b_butter a_butter] = butter(4, 0.2, 'low');
H_butter = freqz(b_butter, a_butter);
 
[b_cheby a_cheby] = cheby1(4, 0.5, 0.2, 'low');
H_cheby = freqz(b_cheby, a_cheby);
 
[b_ellip a_ellip] = ellip(4, 0.5, 40, 0.2, 'low');
H_ellip = freqz(b_ellip, a_ellip);
 
%plot each filter to compare 
figure(11)
norm_freq_axis = [0:1/(512 -1):1];
plot(norm_freq_axis, abs(H_butter))
hold on
plot(norm_freq_axis, abs(H_cheby),'r')
plot(norm_freq_axis, abs(H_ellip),'g')
legend('Butterworth', 'Chebyshev', 'Elliptical')
xlabel('Normalised Frequency');
ylabel('Magnitude')
 
%plot in dB for verification that spec is met
figure(12);
plot(norm_freq_axis, 20*log10(abs(H_butter)))
hold on
plot(norm_freq_axis, 20*log10(abs(H_cheby)),'r')
plot(norm_freq_axis, 20*log10(abs(H_ellip)),'g')
legend('Butterworth', 'Chebyshev', 'Elliptical')
xlabel('Normalised Frequency ');
ylabel('Magnitude (dB)')