function [PAATu,freq]=spec(u,smoothing,fs)

%spectral analysis (gives variance at a particular freqency)
%smoothing is how smooth you want the spectrum
%fs is sampling frequency
%u is data

% number of sampling points per record
num=length(u);
t=1:num;
t=t';

% this sets up where to plot from and to
fr1=1;
to1=num;
i=sqrt(-1);

% setting of the size of the smoothing window
n=smoothing;%num /50;
s=ones(1,2*n+1)/(2*n+1);

% doing the transforms, and removing mean
FAu=fftshift((fft(u(fr1:to1)- mean(u(fr1:to1)))));

% setting up the x axis scale so it reads only to the nyquist
N=length(FAu); 
if(rem(N,2)==1)
  mid=(N-1)/2;
  q=(fs/(2*mid));
  freq=-mid:mid;
  freq=freq*q;
  l=(mid+1):length(FAu);
else 
  mid=(N)/2;
  q=(fs/(2*(mid-1)));
  freq=-mid:(mid-1);
  freq=freq*q;
  l=(mid+1):length(FAu);
end

% calculating the spectra
PAAu=conv(s,FAu.*conj(FAu));            
PAAu=PAAu(n+1:n+N);

%must divide PAA by numer of observations to conserve variance
PAATu=PAAu/num;

%plotting the spectra
% plot(freq(l),PAATu(l),'b')
%   title(['Spectral plot']);
  
  freq=freq(l);
  PAATu=PAATu(l);
  
  