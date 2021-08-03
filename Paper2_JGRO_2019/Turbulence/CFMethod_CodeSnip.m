theta = 30;
%calc along-beam velocities
u1 = -u*sind(theta)-w1*cosd(theta);
u2 = u*sind(theta)-w2*cosd(theta);
u3 = -v*sind(theta)-w1*cosd(theta);
u4 = v*sind(theta)-w2*cosd(theta);

%compute velocity spectra
[Su1,f] = pwelch(u1,[],[],[],fs);
[Su2,~] = pwelch(u2,[],[],[],fs);
[Su3,~] = pwelch(u3,[],[],[],fs);
[Su4,~] = pwelch(u4,[],[],[],fs);

%calculate cospectra
COuw = (Su1-Su2)./(4*cosd(theta)*sind(theta));
COvw = (Su3-Su4)./(4*cosd(theta)*sind(theta));

%Edit: 22/05/18
%flip cospectra u/d so to focus only on the hf
%turbulence now BELOW the waveband (after flipping)
COuw = flipud(COuw);
COvw = flipud(COvw);

%estimate wavenumbers based on U
k = 2*pi*f./U;kc = 2*pi*cf/U;
k = flipud(k); %also flip k (see edit: 22/05/18)
dk = k(2)-k(3);
df = 2*pi*f(3)-2*pi*f(2);

%by definition, the covariance of two signals is the integral of
%the cospectrum (this bit is from Gerbi, 2008 eqn. 6)
uwstar = sum(COuw.*df);vwstar = sum(COvw.*df);
%compute variance-preserving cospectrum
COuwvar = COuw.*k;
COvwvar = COvw.*k; 

%Find the peak of the variance-preserving cospectrum
%above the cutoff frequency
[~,lc] = max(COuwvar);
ku0 = k(lc);
[~,lc] = max(COvwvar);
kv0 = k(lc);
A = (7/3*pi)*sin(3*pi/7);
B = (1/ku0)./(1+((k./ku0).^(7/3)));
COuwstar = uwstar*(A*B);
B = (1/kv0)./(1+((k./kv0).^(7/3)));
COvwstar = vwstar*(A*B);

%compute Ogive curve (CDF) of cospectra
oguw = cumsum(COuw.*dk);
ogvw = cumsum(COvw.*dk);
oguwstar = cumsum(COuwstar.*dk);
ogvwstar = cumsum(COvwstar.*dk);

%turb energy is contained in low f in CDF functions
%lin fit the ogive curves, Edit: 23/05/2018 Steve
%recommends fitting through zero; instead of using
%polyfit (which provides an offset, the y-int), use
%"my_regression" which fits through zero. 

[uw,~,rsq] = my_regression(oguw(k>kc),oguwstar(k>kc),1);
ursq(jj) = rsq;
uw = uw(1);

[vw,~,rsq] = my_regression(ogvw(k>kc),ogvwstar(k>kc),1);
vrsq(jj) = rsq;
vw = vw(1);
