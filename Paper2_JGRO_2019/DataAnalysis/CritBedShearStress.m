%Calculate Tcr, the critical bed shear stress, from Aaron's data!
clear
load('g:\Mekong_W2015\DataAnalysis\Paper3\Sediment\MangroveGrainSizeBen.mat')

NE = Mangroves.MM_098;SW = Mangroves.MM_200;
%Using wet density values from 'The Engineering Toolbox' for clay, silt and
%sand
%wet bulk densities:
sand = 1905;
clay = 1760;
silt = (sand+clay)/2;

%tinoco & coco
p = 1026; %seawater
ps = sand*(SW.sand/100)+clay*(SW.clay/100)+silt*(SW.silt/100);
d50 = 0.001*2^-SW.D50;
sp = ps/p;
g = 9.81;
nu = 1.05E-6;
Ds = (((g*(sp-1))/nu^2)^(1/3))*d50;
theta = (0.3/(1+1.2*Ds))+0.055*(1-exp(-0.020*Ds));
Tcr = theta*g*(ps-p)*d50;


