function out = compute_spectra(time,data,win,step,w,fs);
%Compute auto spectra and cross spectra of variables for H, Bd, and X, Bd
%analysis
% 
% Inputs: time - vector of times
%         data - structure containing data to compute spectra from
%          win - window length for spectral analysis
%         step - step interval between windowed segments
%            w - spectra window length
%           fs - sampling frequency in Hz
%
% Outputs: out - structure containing spectral/cross spectral analysis
%
% Paper 3: Sediment Motion in Mangroves
%
% This is Version 1.0 of this script
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%inputs
fn = fieldnames(data);
if ~any(strcmp(fn,'h')) && ~any(strcmp(fn,'bd')) && ~any(strcmp(fn,'x'))
    disp('Cannot run without depth, bottom trace and x-shore vels')
    return;
end

avt = step*fs; %# samples/step
nwin = win*fs; %# samples/window
for i = 1:length(fn);
    nsamp = length(data.(fn{i}));
    ind = [1 avt:avt:nsamp];
    for j = 1:length(ind)
        if abs(nsamp-ind(j)) < nwin  %skip the last few indexes approaching the end of the t-s
            continue
        else
            idx = ind(j):ind(j)+nwin-1;
        end
        t = time(idx(1));
        var = detrend(data.(fn{i})(idx));
        N = length(var);
        M = w;
%         nu = 2*floor(N/M); %DOF
        sw = hamming(M);
        [psd,f,pxxc]=pwelch(var,sw,[],N,fs,'confidencelevel',0.95);
%         alfa = 1 - 0.95;
%         c = chi2inv([1-alfa/2 alfa/2],nu);
%         c = nu./c; %spectra confidence intervals
        pxxc = mean(pxxc);
        %save to structure
        if j==1,out.(fn{i}).f = f;end
        out.(fn{i}).psd(j,:) = psd;
        out.(fn{i}).ci(j,:) = pxxc;
        out.(fn{i}).time(j,:) = t;
    end
end

%compute power spectra for h-bd and x-bd
for i = 1:2
    if i == 1
        var1 = data.h;
        fn = 'hbd';
    else
        var1 = data.x;
        fn = 'xbd';
    end
    var2 = data.bd;
    
    nsamp = length(var1);
    ind = [1 avt:avt:nsamp];
    for j = 1:length(ind)
        if abs(nsamp-ind(j)) < nwin  %skip the last few indexes approaching the end of the t-s
            continue
        else
            idx = ind(j):ind(j)+nwin-1;
        end
        v1 = detrend(var1(idx)); %h or x
        v2 = detrend(var2(idx)); %bd
        
        [s12,f] = cpsd(v2,v1,sw,[],N,fs);
        [s1,~] = pwelch(v1,sw,[],N,fs);
        [s2,~] = pwelch(v2,sw,[],N,fs);
        coh = abs(s12).^2 ./ (s1.*s2);
        phase= atan2(-imag(s12),real(s12))*180/pi;
        nu = 2*floor(N/M); %DOF
        ci = cohere_signif_level(nu);
        if j==1,out.(fn).f = f;end
        out.(fn).cpsd(j,:) = s12;
        out.(fn).coh(j,:) = coh;
        out.(fn).phase(j,:) = phase;
        out.(fn).ci(j,:) = ci;
    end
end
