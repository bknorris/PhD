
function tce = critstresswu(psand,psilt,pclay,d50,type)
%Calculate estimated critical shear stress for mud + sand mixtures (Wu et
%al. 2017)
%
% Inputs: psand, psilt, pclay: percentage sand, silt and clay (up to 1)
%         d50: median grain size in m.
%         type: either 'cohesive' or 'noncohesive'. Cohesive sediments are
%         >10% mud, and non-cohesives <10% mud.
% Outputs: tce, the critical erosion threshold of a mud+sand mixture based
%          off two datasets: Smith (2015) for cohesives, and 
%          Panagiotopoulos et al. 1997 for non-cohesives. 
%
% See Wu et al. 2012 for details. Code by BKN, Waikato 2017
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
pmud = psilt+pclay;
rhod = 1100; %dry density of mud+sand mixture (our data)
betav = 1.2; %beta_v set to 1.2 (Wu et al. 2017)
tcrn = 0.15; %N/m^2, from Shields (1936) and Soulsby (1997)
if strcmp(type,'cohesive')
    %Smith 2015 data
    rhos = 2600; %dry density of sand (Ye et al. 2011)
    rhom = 2575; %dry density of mud (Ye et al. 2011)
    ds = 0.353/1000; %sand diameter (Smith, 2015)
    dm = median([0.0038 0.005]./1000); %mud diameter (Smith, 2015)
    Bmax = 0.65; %Bmax from Wu et al. (2017)
    tce = [0.544 1.403 1.403 1.983 1.395];
    rc = [0.296 0.253 0.253 0.364 0.339];
    alpha = [0.08 0.09 0.1 0.17 0.08];
    %Calculations
    Rn = (pmud*ds^3)/(psand*dm^3); %eq. 23
    P = (pmud/dm)/((psand/ds)+(pmud/dm)); %eq. 24
    n = 0.1124*(ds/dm); %eq. 25
    Nc = pi/(asin(dm/(ds+dm))^2)*cosd(30); %eq. 26
    B = min([(P*Rn)/(n*Nc),Bmax]); %eq. 22
    rhodm = pmud/((1/rhod)-((psand/rhos)*B)+((psand/rhos)*(1-B))); %eq. 21
    phim = 1-(rhodm/rhos);
    r = (1-phim)/phim;
    tildece = tce.*(r./rc).^1.7;
    tcel = tcrn+1.25*(tce-tcrn)*min([pmud,0.05]);
    tcem = tcel+(tildece-tcel).*exp(-alpha.*(psand/pmud)^betav);
elseif strcmp(type,'noncohesive')
    %Panagiotopoulos et al. 1997 data
    dm = 0.0032/1000; %mud diameter
    ds = [0.1525 0.215]; %sand diameter
    rhos = [1481 1538]; %dry density of sand (estimated, Wu et al. 2017)
    rhom = 892; %dry mud density, estimated using formula from Wu & Wang, 2006
    Bmax = 0.65;
    tce = 0.144;
    rc = [0.232 0.241];
    alpha = [0.275 0.24];
    
    
    
    
end