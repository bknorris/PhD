%%%%%%%%%%%%%%%%%
% DFT exaplained
%
% SMJ
%
% 3/21/2014
%
%%%%%%%%%%%%%%%%%
close all; clear
n=0:2^5-1;
N = length(n);
% single cosine
y = 3*cos(2*pi*n/N);
% two waves
%y = 4*sin(2*pi*1*n/N) + 2*cos(2*pi*3*n/N);
% wave with a phase shift: 2.6525sin(2pi(2)n/N) + 1.4cos(2pi(2)n/N) 
%y = 3*cos(2*pi*2*n/N - 1.0854);
% noniteger number of cycles
%y = 5*cos(2*pi*2.3*n/N);
%y = [y, zeros(1,1000)];
N = length(y);
n = 0:N-1;
%% plot a test function
orient tall
figure
subplot(10,2,1:2)
p1=plot(n,y,'.-');
set(p1,'MarkerFaceColor','r');
set(p1,'MarkerEdgeColor','r');
grid on

%% compute some basis functions
basis_cos0 = cos(2*pi*0*n/N);
basis_sin1 = sin(2*pi*1*n/N);
basis_cos1 = cos(2*pi*1*n/N);
basis_sin2 = sin(2*pi*2*n/N);
basis_cos2 = cos(2*pi*2*n/N);
basis_sin3 = sin(2*pi*3*n/N);
basis_cos3 = cos(2*pi*3*n/N);

%% plot basis functions
plot_basis_functions(n,basis_cos0,basis_sin1,basis_cos1,basis_sin2,basis_cos2,basis_sin3,basis_cos3)

%% take convolution of test function * basis function
conv_cos0 = y .* basis_cos0;
conv_sin1 = y .* basis_sin1;
conv_cos1 = y .* basis_cos1;
conv_sin2 = y .* basis_sin2;
conv_cos2 = y .* basis_cos2;
conv_sin3 = y .* basis_sin3;
conv_cos3 = y .* basis_cos3;

%% plot convolution
plot_convolution(n,conv_cos0,conv_sin1,conv_cos1,conv_sin2,conv_cos2,conv_sin3,conv_cos3)

%% Amplitudes:  divide summed convolution by N/2