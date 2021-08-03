function z0 = baptistroughness(H,Cd,hc,N,d,phi)
% Calculate z0 from Baptist et al. (2007) based on the vegetation
% characteristics, water depth and cp (the turbulence closure coefficient). 
% see Baptist et al. (2007) for details.
% 
% Inputs: H - the water depth (m)
%         Cd - the canopy drag coefficient
%         hc - the canopy height (m)
%         N - number of stems (m^-2)
%         d - mean stem diameter (m)
%         phi - volume fraction occupied by vegetation
%
% Outputs: z0 - bed roughness length scale (for shear stress calculations)
%
% Code by BKN, Waikato 2017
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%define terms
cl = 1;
l = cl*((1-phi)/N)^(1/2); %eq. 25
cp = (1/20)*((H-hc)/l); %cp determined by a fit of Nepf & Vivoni 2000 data. 
%Needs to be calibrated with other data!
L = sqrt((cp*l)/(Cd*N*d)); %eq. 31
D = hc-L*(1-exp((-1*hc/L))); %eq. 47
z0 = (hc-D)*exp(-0.4*sqrt((2*L/cp*l)*(1+(L/(H-hc))))); %eq. 49