function plot_convolution(n,conv_cos0,conv_sin1,conv_cos1,conv_sin2,conv_cos2,conv_sin3,conv_cos3)

%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here
subplot(8,2,4)
plot(n,conv_cos0,'.-','MarkerFaceColor','r','MarkerEdgeColor','r');
hold on
plot(n,zeros(length(n),1),'k')
title(['Sum = ',num2str(sum(conv_cos0))])

subplot(8,2,6)
plot(n,conv_sin1,'.-','MarkerFaceColor','r','MarkerEdgeColor','r');
hold on
plot(n,zeros(length(n),1),'k')
title(['Sum = ',num2str(sum(conv_sin1))])

subplot(8,2,8)
plot(n,conv_cos1,'.-','MarkerFaceColor','r','MarkerEdgeColor','r');
hold on
plot(n,zeros(length(n),1),'k')
title(['Sum = ',num2str(sum(conv_cos1))])

subplot(8,2,10)
plot(n,conv_sin2,'.-','MarkerFaceColor','r','MarkerEdgeColor','r');
hold on
plot(n,zeros(length(n),1),'k')
title(['Sum = ',num2str(sum(conv_sin2))])

subplot(8,2,12)
plot(n,conv_cos2,'.-','MarkerFaceColor','r','MarkerEdgeColor','r');
hold on
plot(n,zeros(length(n),1),'k')
title(['Sum = ',num2str(sum(conv_cos2))])

subplot(8,2,14)
plot(n,conv_sin3,'.-','MarkerFaceColor','r','MarkerEdgeColor','r');
hold on
plot(n,zeros(length(n),1),'k')
title(['Sum = ',num2str(sum(conv_sin3))])

subplot(8,2,16)
plot(n,conv_cos3,'.-','MarkerFaceColor','r','MarkerEdgeColor','r');
hold on
plot(n,zeros(length(n),1),'k')
title(['Sum = ',num2str(sum(conv_cos3))])

end

