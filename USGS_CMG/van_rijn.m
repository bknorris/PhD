function [LAMvr ETAvr] = van_rijn(T,d_50,Ub_rms)

g=9.81 *100; %m/s2 -->cm/s
density_water = 1.028; %g/cm3 (seawater)
density_sed = 2.650; %g/cm3
s=density_sed/density_water;
d=d_50.*1e-4; % microns --> cm
%Ub_sig = (1.42 * sqrt(2) * Ub_rms);

mobility=Ub_rms.^2/(g*s*d);

omega = (2*pi)./T; %1/s

A=(Ub_rms)./omega * (100); %m -->cm

Ia=find(mobility <= 10);
    ETA_over_a(Ia) = 0.22;
    ETA(Ia)=ETA_over_a(Ia).*A(Ia);
    ETA_over_LAM(Ia) = 0.18;
    LAM(Ia)=ETA(Ia)/ETA_over_LAM(Ia);
    
Ib=find(mobility <250 & mobility >10);
    ETA_over_a(Ib) = 2.8e-13 .* (250- mobility(Ib)).^5;
    ETA(Ib)=ETA_over_a(Ib).*A(Ib);
    ETA_over_LAM(Ib) = 2e-7.*(250-mobility(Ib)).^2.5;
    LAM(Ib)=ETA(Ib)/ETA_over_LAM(Ib);
    
Ic=find(mobility >= 250);
    ETA_over_a(Ic) = 0;
    ETAvr(Ic)=ETA_over_a(Ic).*A(Ic);
    %ETA_over_LAM(Ic) = 0;
    LAMvr(Ic)=0;
    
Ic=find(isnan(mobility));
    ETAvr(Ic)=NaN;
    %ETA_over_LAM(Ic) = NaN;
    LAMvr(Ic)=NaN;


    