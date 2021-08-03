function k=om2k(omega,h);

%function k=om2k(omega,h);
%
% matlab function to invert surface gravity wave dispersion relation, 
% i.e. to calculate radian wavenumber k (m^{-1}) from 
% radian frequency omega (s^{-1}) and depth h (m).

warning off MATLAB:fzero:UndeterminedSyntax
g=9.81;
k=nan*ones(size(omega));
for ii=1:length(omega)
    if isnan(omega(ii)*h)
        k(ii)=nan;
    else
        if omega(ii)==0
            k(ii)=0;
        else
            k0=omega(ii)^2/g;
            k0=k0/sqrt(tanh(k0*h));
            k(ii)=fzero('dispersion',k0,[],abs(omega(ii)),h);
        end
    end
end