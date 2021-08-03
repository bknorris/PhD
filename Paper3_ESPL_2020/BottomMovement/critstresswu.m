function tcem = critstresswu(psand,pmud,d50,type)
%Calculate estimated critical shear stress for mud + sand mixtures
%
% Inputs: psand, psilt, pclay: percentage sand and mud (equal to 1)
%         d50: median grain size in m.
%         type: either 'muddy' or 'sandy'. Muddy is defined as sediments
%         with < 15% sand, and Sandy as sediments with > 15% sand.
% Outputs: tce, the critical erosion threshold of a mud+sand mixture based
%          off two datasets: Smith (2015) for cohesives, and
%          Panagiotopoulos et al. 1997 for non-cohesives.
%
% See Wu et al. 2012 for details. Code by BKN, Waikato 2017
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
betav = 1.2; %beta_v set to 1.2 (Wu et al. 2017)
tcrn = 0.12; %N/m^2, from Shields (1936) and Soulsby (1997)
pmapos = 0.13+(0.21/(d50*1000+0.002)^0.21); %eq 4 from Wu & Wang, 2006
rhowm = 1840; %wet density of mud (Eng. Toolbox)
rhows = 1905; %wet density of sand (Eng. Toolbox)
rhods = 1555; %dry density of sand (Eng. Toolbox)
if strcmp(type,'muddy')
    %Smith 2015 data
    rhod = (1-pmapos)*rhowm; %dry mixture density based off of d50.
    dm = d50; %mud diameter [m]
    ds = 0.000353; %sand diameter [m] (Smith, 2015)
    Bmax = 0.73; %Bmax from Wu et al. (2017) - mean value
    tce = [0.544 1.403 1.403 1.983 1.395];
    rc = [0.296 0.253 0.253 0.364 0.339];
elseif strcmp(type,'sandy')
    %Panagiotopoulos et al. 1997 data
    rhod = (1-pmapos)*rhows; %dry mixture density based off of d50.
    dm = 0.0000127; %mud diameter [m] (Smith, 2015)
    ds = d50; %sand diameter [m]
    Bmax = 0.73; %Bmax from Wu et al. (2017) - mean value
    tce = [0.144 0.144];
    rc = [0.232 0.241];
end
%Calculations
Rn = (pmud*ds^3)/(psand*dm^3); %eq. 23
P = (pmud/dm)/((psand/ds)+(pmud/dm)); %eq. 24
n = 0.1124*(ds/dm); %eq. 25
Nc = pi/(asin(dm/(ds+dm))^2)*cosd(30); %eq. 26
B = min([(P*Rn)/(n*Nc),Bmax]); %eq. 22
rhos = (psand*rhows)+(pmud*rhowm);
rhodm = (rhod*rhos*rhods*pmud)/(rhod*psand*(rhos*(B-1)-rhods*B)+rhos*rhods); %eq. 21, solve for rhodm
phim = 1-(rhodm/rhos);
r = (1-phim)/phim;
tildece = tce.*((r./rc).^1.7); %eq. 20
if pmud <= 0.05
    tcel = tcrn+1.25*(tce-tcrn)*min([pmud,0.05]); %eq. 17
else
    tcel = tce;
end
alpha = 0.42*exp(-3.38*ds); %eq. 27
tcem = min(tcel+(tildece-tcel).*exp(-alpha.*(psand/pmud)^betav)); %eq. 19
end