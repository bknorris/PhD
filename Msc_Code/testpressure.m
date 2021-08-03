load S80Pburst.mat
close all
P=detrend(S80P.p);

F=fft(P);
fn=1; %nyquist frequency
f=0:599; f=f./599;
f=[f f(600:-1:1)];
%

subplot(221)
plot(f,abs(F))
d=1.0892;
z=-d;

[kk]=dispersion(f(1:600),d);
k=2*pi*kk;
k=[k k(600:-1:1)];


attn=cosh(k.*(d+z))./cosh(k.*d);


subplot(222)


plot(f,attn)

FF=F./attn;

j=attn<0.3;
FF(j)=0;

subplot(223)
plot(f,abs(FF),'r')

PP=real(ifft(FF));

subplot(224)
hold off
plot(P)
hold on
plot(PP,'r')
