function plot_basis_functions(n,basis_cos0,basis_sin1,basis_cos1,basis_sin2,basis_cos2,basis_sin3,basis_cos3)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
subplot(8,2,3)
plot(n,basis_cos0,'.-','MarkerFaceColor','r','MarkerEdgeColor','r');
title('cos(2\pi(0)n/N')

subplot(8,2,5)
plot(n,basis_sin1,'.-','MarkerFaceColor','r','MarkerEdgeColor','r');
title('sin(2\pi(1)n/N')

subplot(8,2,7)
plot(n,basis_cos1,'.-','MarkerFaceColor','r','MarkerEdgeColor','r');
title('cos(2\pi(1)n/N')

subplot(8,2,9)
plot(n,basis_sin2,'.-','MarkerFaceColor','r','MarkerEdgeColor','r');
title('sin(2\pi(2)n/N')

subplot(8,2,11)
plot(n,basis_cos2,'.-','MarkerFaceColor','r','MarkerEdgeColor','r');
title('cos(2\pi(2)n/N')

subplot(8,2,13)
plot(n,basis_sin3,'.-','MarkerFaceColor','r','MarkerEdgeColor','r');
title('sin(2\pi(3)n/N')

subplot(8,2,15)
plot(n,basis_cos3,'.-','MarkerFaceColor','r','MarkerEdgeColor','r');
title('cos(2\pi(3)n/N')

end

