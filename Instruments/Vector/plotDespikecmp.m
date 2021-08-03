figure
ax1 = subplot(311);
plot(B1,'color','k')
hold on
plot(V1,'r')
legend('No Despike','With Despike')
title('V1')
xlabel('N samples')
ylabel('m/s')
ax2 = subplot(312);
plot(B3,'k')
hold on
plot(V2,'r')
legend('No Despike','With Despike')
title('V2')
xlabel('N samples')
ylabel('m/s')
ax3 = subplot(313);
plot(B3,'k')
hold on
plot(V3,'r')
legend('No Despike','With Despike')
title('V3')
xlabel('N samples')
ylabel('m/s')
linkaxes([ax1 ax2 ax3],'x')
% axis([0.3 2.1E6 -1.5 1.5])