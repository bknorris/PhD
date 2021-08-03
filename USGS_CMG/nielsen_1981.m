function [LAMn_1981,ETAn_1981]=nielsen_1981(T,d_50,ubr)
%Nielsen (1981) ripple prediction
%need: wave period, median grain diameter, Ub_sig, sediment density
%Created By Tim Nelson
%University of South Carolina, Department of Geological Sciences
%adjusted for units by Jordan Landers (WHOI-SSF09) 7/7/2009

g=9.81; %m/s2
density_water = 1028; %kg/m3 (seawater)
density_sed = 2650; %kg/m3
s=density_sed/density_water;
Ub_sig=sqrt(2).*ubr;
omega = (2*pi)./T; %1/s
d=d_50.*1e-6;

%Ub_sig= 1.42 * sqrt(2) * Ub_rms;
a=(Ub_sig)./omega; %m
                                 
Mobility_number = ((a.*omega).^2)./((s-1).*g.*d);

% Friction_Factor = exp(5.213 .* (((2.5 * d) ./ a).^0.194) - 5.977);

%Shields_Parameter = 0.5 .* Friction_Factor .* Mobility_number;

lambda_over_a = exp((693 - 0.37 .* (log(Mobility_number)).^8)./(1000 + 0.75 .* (log(Mobility_number)).^7));

lambda = lambda_over_a .* a;

LAMn_1981 = lambda; %m

%steepness = (0.342 - 0.34 .* (Shields_Parameter).^0.25)/100;


j=find(Mobility_number >= 10);
k=find(Mobility_number < 10);
l=find(isnan(Mobility_number));


height_over_a(j) = 21 .* Mobility_number(j).^-1.85;  %only valid for quartz sand
height_over_a(k) = 0.275 - 0.022 .* sqrt(Mobility_number(k));
height_over_a(l) = NaN;

ETAn_1981 = (height_over_a(:) .* a(:)); %in m