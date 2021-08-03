%Patch wake length scale for the HTA experiment, given Chen's (2012)
%empirical relationship for Lw (wake length) and D (patch width).

% Cd = 2; %from Henderson et al. (2016)
a = 4; %from veg reconstructions (HTA day 1)
D = 0.21; %measured in SelectedPneumatophoreStatistics
Cd = 2; %from Fig. 5 of Tanino & Nepf
p = [-10 0 10];
q = [-0.2 0 0.2];
[P,Q] = meshgrid(p,q);
c=cat(2,P',Q');
d=reshape(c,[],2); %all combinations of p and q
Lw_D = zeros(length(d),1);
for i = 1:length(d)
    n = 25+d(i,1);
    m = -0.9+d(i,2);
    Lw_D(i) = 1.2+n*(Cd*a*D)^(m);
end
L1 = (10*((2-0.25*Cd*a*D)/Cd*a*D))*D;
disp(['CdaD is: ' sprintf('%0.3f',Cd*a*D)])
disp(['Distance to onset of steady wake ' sprintf('%0.3f',L1) ' m'])
Lw = min(Lw_D)*D;
disp(['Minimum patch wake length: ' sprintf('%0.3f',Lw) ' m'])
Lw = mean(Lw_D)*D;
disp(['Mean patch wake length: ' sprintf('%0.3f',Lw) ' m'])
Lw = max(Lw_D)*D;
disp(['Maximum patch wake length: ' sprintf('%0.3f',Lw) ' m'])