%Learn wavelet coherence analysis

%Example from the File Exchange
d1=load('jao.txt'); 
d2=load('jbaltic.txt'); 
wtc(d1,d2)

%Trying to build an example from Mathworks
%Sample rate is 10 and 50Hz for x and y, respectively

tx = 0:0.01:2;
x = cos(2*pi*10*tx).*(tx>=0.5 & tx<1.1)+ ...
cos(2*pi*50*tx).*(tx>= 0.2 & tx< 1.4)+0.25*randn(size(tx));
y = sin(2*pi*10*tx).*(tx>=0.6 & tx<1.2)+...
sin(2*pi*50*tx).*(tx>= 0.4 & tx<1.6)+ 0.35*randn(size(tx));

figure
plot(tx,x,'b'),hold on
plot(tx,y,'r')

%try some things; the following parameters are from Torrence & Compo, 1998
fs = 20;
dt = 1/fs;
s0 = 2*dt;
Dj = 1/12;
N = length(x);
J1 = (Dj^-1)*log2(N*dt/s0);
[Rsq,period,scale,coi,sig95,Wxy,t,dt]=wtc(x,y,'S0',10/60);
plotwtc(Rsq,period,coi,sig95,t,Wxy,dt)