figure('units','normalized','outerposition',[0 0 1 1])
a = subplot(311);
plot(ADV.datetime,ADV.Pres,'k')
datetickzoom('x','dd HH:MM:SS','keepticks')
ylabel('\bf\itdBar')
title('\bfADV Pressure Signal')
b = subplot(312);
plot(ADV.datetime,ADV.V1,'b')
hold on
plot(ADV.datetime,ADV.V2,'r')
hold on
plot(ADV.datetime,ADV.V3,'g')
datetickzoom('x','dd HH:MM:SS','keepticks')
ylabel('\bf\itm/s')
title('\bfUnmodified Velocities X Y Z')
legend('\bf\itV1','\bf\itV2','\bf\itV3')
c = subplot(313);
plot(ADV.datetime,ADV.U,'b')
hold on
plot(ADV.datetime,ADV.V,'r')
hold on
plot(ADV.datetime,ADV.W,'g')
datetickzoom('x','dd HH:MM:SS','keepticks')
ylabel('\bf\itm/s')
title('\bfRotated Velocities U V W')
legend('\bf\itU','\bf\itV','\bf\itW')
set([a b c],'XGrid','on')
xlabel('\bf\itdd HH:MM:SS')
linkaxes([a,b,c],'x')