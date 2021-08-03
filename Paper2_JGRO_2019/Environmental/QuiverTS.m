N = 56;
Up = 360-[sw.day1.udir sw.day2.udir sw.day3.udir];
USp = [sw.day1.uspd sw.day2.uspd sw.day3.uspd];
tt = date2doy([sw.day1.time sw.day2.time sw.day3.time])';
[u,v] = cmgspd2uv(USp,Up);
ys = zeros(length(tt),1);
subplot(211)
quiver(tt,ys,u'.*0.1,v'.*0.1,'AutoScale','off')

Lw = 360-[sw.day1.ldir sw.day2.ldir sw.day3.ldir];
LSp = [sw.day1.lspd sw.day2.lspd sw.day3.lspd];
[u,v] = cmgspd2uv(LSp,Lw);
ys = zeros(length(tt),1);
subplot(212)
quiver(tt,ys,u'.*0.1,v'.*0.1)
