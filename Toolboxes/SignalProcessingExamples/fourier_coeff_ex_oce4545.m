%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% OCE 5586
%
% Fourier Coefficients example
%
% SMJ
% 2/27/2014
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
close all; clear;

% actual data
n_sst = 1:24;
sst = [7.6 7.4 8.2 9.2 10.2 11.5 12.4 13.4 13.7 11.8 10.1 9.0 ...
    8.9 9.5 10.6 11.4 12.9 12.7 13.9 14.2 13.5 11.4 10.9 8.1];

figure
plot(n_sst,sst)
grid
xlabel('n')
ylabel('SST (deg C)')


% get Fourier coeffient, A0
%
N = length(n_sst);

A0 = 0;
for n = 1:N
    A0 = A0 + sst(n); 
end
A0 = 2/N*A0;
B0 = 0;

% get Fourier coeffients, Ap and Bp
%
A = zeros(N/2,1);
B = zeros(N/2,1);

for p = 1:N/2
    for n = 1:N
        A(p) = A(p) + sst(n) * cos(2*pi*p*n/N);
        B(p) = B(p) + sst(n) * sin(2*pi*p*n/N);
    end
end
A = 2/N*A;
B = 2/N*B;

% reconstruction
for max_index=1:10
sst_reconstr = zeros(N,1);
for n = 1:N
    for p = 1:max_index
        sst_reconstr(n) = sst_reconstr(n) + A(p)*cos(2*pi*p*n/N) ...
            + B(p)*sin(2*pi*p*n/N);
    end
end
sst_reconstr = A0/2 + sst_reconstr;

% figure
hold on
plot(n_sst,sst_reconstr,'r')
legend('Actual','Fourier fit')
pause
end
        
% prediction
% sst_reconstr2 = zeros(48,1);
% for n = 1:48
%     for p = 1:3
%         sst_reconstr2(n) = sst_reconstr2(n) + A(p)*cos(2*pi*p*n/N) ...
%             + B(p)*sin(2*pi*p*n/N);
%     end
% end
% sst_reconstr2 = A0/2 + sst_reconstr2; 
% 
% plot(1:48,sst_reconstr2,'g')