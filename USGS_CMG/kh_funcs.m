function kh = kh_funcs(w,h,method)
% KH_FUNCS - Various functions for solving wave-dispersion relation
% function kh = kh_funcs(w,h,method)
%
% Hardwired for MKS units with g = 9.80665 m/s^2
%
% Input:
%  w = radian wave frequency = 2*pi/T with T = period (s)
%  h = water depth (m)
%  method (optional)
%     'iterate'  - Newton-Raphson iteration to spec. accuracy in kh
%     'iterate2' - Newton-Raphson iteration using Soulsby, 2006 method
%     'gilbert'  - G. Gilbert method cited in Soulsby, 2006, Eqn. 17-18
%     'fenton'   - Fenton and McKee (1989) as quoted by Fenton (1990)
%                   The Sea, vol. 9A, and Dalrymple's web page
%     'carvalho1'- Method derived by R. Carvalho using Gene expression
%                   programming (Eqn. 1 in his attachment to Coastal list)
%     'carvalho2'- Method derived by R. Carvalho using Gene expression
%                   programming (Eqn. 2 in his attachment to Coastal list)
%     'carvalho3'- Method derived by R. Carvalho using Gene expression
%                   programming (Average of his Eqn. 1 + 2)
%     'hunt'     - (Default method) Hunt (1979) cited in Dean and
%                   Dalrymple (1991) p. 72
% Returns:
%  kh = wavenumber * depth with k = 2*pi/L and L = wavelength (m)
%
% Note: We found that the if/else construct and strcmpi calls slowed these
% routines down by almost tenfold, and the routines lower in the if/else
% cascade were slower than the ones near the top. 
% We recommend implementing either Hunt or iterate2 as
% a stand-alone function.

% csherwood@usgs.gov
% Last revised Sept. 10, 2006
g = 9.80665;

if(exist('method','var')~=1), method='default'; end

if( (strcmpi(method,'hunt'))||(strcmpi(method,'default')))
   x = (w.^2) .* (h ./g);
   kh2 = 1.0 + x.*(0.6666666666+x.*(0.3555555555+x.*(0.1608465608+...
      x.*(0.0632098765+x.*(0.0217540484+x.*0.0065407983)))));
   kh = sqrt(x.^2+x./kh2);   
elseif(strcmpi(method,'iterate'))
   % Solve:  0 = kh*tanh(kh)-w**2*h/g using Newton's method
   jmax=20;
   xacc=0.000001; % normally 0.0001
   w2h = w*w*h/g;
   %...BOUNDS ARE 0 AND SLIGHTLY MORE THAN W2
   x1=0.;
   x2=w2h+1.;
   %...BEST INITIAL GUESS IS W2
   kh=w2h;
   for j=1:jmax,
      %       ...COMPUTE F(X)
      tx=tanh(kh);
      f=kh*tx-w2h;
      %       ...COMPUTE F'(X)
      df=kh-kh*tx*tx + tx;
      dx=f/df;
      kh=kh-dx;
      if((x1-kh)*(kh-x2)<0),
         error('kh_funcs solver out of brackets');
      end;
      if(abs(dx)<xacc),break,end;
   end
   if(j==jmax),error('kh_funcs exceeded maximum iterations'),end
elseif( strcmpi(method,'iterate2'))
   x = w.^2*h./g;
   y = sqrt(x) .* (x<1) + ...
      x.* (x>=1);
   %for k=1:3,
   t = tanh( y );
   y = y-( (y.*t -x)./(t+y.*(1-t.^2)));
   t = tanh( y );
   y = y-( (y.*t -x)./(t+y.*(1-t.^2)));
   t = tanh( y );
   y = y-( (y.*t -x)./(t+y.*(1-t.^2)));
   %end
   kh=y;
elseif( strcmpi(method,'gilbert'))
   x = w.^2*h./g;
   y = (sqrt(x).*(1+0.2*x)) .* (x<=1) + ...
      x.*(1+ 0.2.*exp(2-2*x)) .* (x>1);
   kh = y;
elseif(strcmpi(method,'fenton'))
   T = 2*pi ./w;
   Ldeep = g*(T.^2)./(2*pi);
   L = Ldeep.*(tanh( (w.^2 .* h / g).^(3/4) )).^(2/3);
   k = 2*pi ./L;
   kh = k.*h;
elseif( strcmpi(method,'carvalho1'))
   T = 2*pi ./w;
   Ldeep = g*(T.^2)./(2*pi);
   k0h = h.*(2*pi./Ldeep);
   L = Ldeep.*tanh( k0h ./ ((tanh(k0h)).^(1/4) .* ...
      sqrt(tanh(sqrt(sinh(k0h))))));
   k = 2*pi ./L;
   kh = k.*h;
elseif( strcmpi(method,'carvalho2'))
   T = 2*pi ./w;
   Ldeep = g*(T.^2)./(2*pi);
   k0h = h.*(2*pi./Ldeep);
   L = Ldeep.*tanh( k0h ./ ...
      tanh( k0h./ ...
      tanh( k0h./ ...
      sinh( tanh(sqrt(k0h))))));
   k = 2*pi ./L;
   kh = k.*h;
elseif( strcmpi(method,'carvalho3'))
   T = 2*pi ./w;
   Ldeep = g*(T.^2)./(2*pi);
   k0h = h.*(2*pi./Ldeep);
   L = Ldeep.* (0.5*tanh( k0h ./ ((tanh(k0h)).^(1/4) .* ...
      sqrt(tanh(sqrt(sinh(k0h)))))) +...
      0.5*tanh( k0h ./ ...
      tanh( k0h./ ...
      tanh( k0h./ ...
      sinh( tanh(sqrt(k0h)))))));
   k = 2*pi ./L;
   kh = k.*h;
else
   fprintf(1,'Method = %s\n',method);
   error('Bad value of method in kh_funcs');
end