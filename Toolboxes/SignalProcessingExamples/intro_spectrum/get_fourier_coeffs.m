function [A0,A,B,N] = get_fourier_coeffs(data_array,plotit);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% OCE 5586
%
% Get Fourier Coefficients
%
% SMJ
% 3/13/2014
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
data = data_array;

%% get length of data
N = length(data);

%% get Fourier coeffient, A0

A0 = 0;
for n = 1:N
    A0 = A0 + data(n); 
end
A0 = 2/N*A0;
B0 = 0;

%% get Fourier coeffients, Ap and Bp
A = zeros(N/2,1);
B = zeros(N/2,1);

for p = 1:N/2
    for n = 1:N
        A(p) = A(p) + data(n) * cos(2*pi*p*n/N);
        B(p) = B(p) + data(n) * sin(2*pi*p*n/N);
    end
end
A = 2/N*A;
B = 2/N*B;

%% plot coefficents
if (plotit==1)
    figure
    subplot(3,1,1)
    bar(A)
    grid
    ylabel('A Coeff')
    subplot(3,1,2)
    bar(B,'g')
    grid
    ylabel('B Coeff')
    subplot(3,1,3)
    bar(A.^2+B.^2,'r')
    grid
    ylabel('A^2 + B^2')
    xlabel('Index, p')
end

% % reconstruction
% data_reconstr = zeros(N,1);
% for n = 1:N
%     for p = 1:N/2
%         data_reconstr(n) = data_reconstr(n) + A(p)*cos(2*pi*p*n/N) ...
%             + B(p)*sin(2*pi*p*n/N);
%     end
% end
% data_reconstr = A0/2 + data_reconstr;
% 
% % figure
% hold on
% plot(time,data_reconstr,'r')
% legend('Actual','Fourier fit')
% 
% return
% 
% % prediction
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