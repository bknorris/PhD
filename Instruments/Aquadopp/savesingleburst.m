%Get a single burst
load('5116_ver3.mat')
spb = aqdp.metadata.spb;
burstl = (length(aqdp.burst)/spb);
burstc = ceil(burstl);
corav = (aqdp.cor1+aqdp.cor2+aqdp.cor3)/3;
backav = (aqdp.beam1+aqdp.beam2+aqdp.beam3)/3;
disp('Pick a Burst')
figure(7)
imagesc(aqdp.datenum,aqdp.rangebins,aqdp.cor1')
datetick('x','keepticks')
ylabel('%')
title('Aquadopp Correlations')
axis xy
hold on
q = ginput(2);
int = find(aqdp.datenum >= q(1,1) & aqdp.datenum <= q(2,1));
s = ceil(int(1)/spb);f = ceil(int(end)/spb);

ix1 = spb*s-(spb-1);
ix2 = spb*f; %this value must always be n+4096 > than ix1
ind = ix1:ix2;
burst = struct();
burst.u = aqdp.u(ind,:);
burst.v = aqdp.v(ind,:);
burst.w = aqdp.w(ind,:);
burst.corav = corav(ind,:);
burst.backav = backav(ind,:);
burst.cor1 = aqdp.cor1(ind,:);
burst.cor2 = aqdp.cor2(ind,:);
burst.cor3 = aqdp.cor3(ind,:);
burst.beam1 = aqdp.beam1(ind,:);
burst.beam2 = aqdp.beam2(ind,:);
burst.beam3 = aqdp.beam3(ind,:);
burst.datenum = aqdp.datenum(ind,:);
burst.yearday = aqdp.yearday(ind,:);
burst.rangebins = aqdp.rangebins;
burst.pressure = aqdp.pressure(ind,:);
burst.nBursts = s:f;
save('5116_single.mat','burst')
close