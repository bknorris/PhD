function out = compute_spectra(data,win,step,w,fs);
%Compute auto spectra and cross spectra of variables for H, Bd, and X, Bd
%analysis
% 
% Inputs: data - structure containing data to compute spectra from
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
        
        var = detrend(data.(fn{i})(idx));
        N = length(var);
        M = w;
        nu = 2*floor(N/M); %DOF
        subid = [1 M/2:M/2:N]; %50% overlap
        sl = length(subid)-2;
        sw = hamming(M+1); %hamming window
        fd = zeros(M+1,sl); %preallocate
        for n=1:sl
            sid = subid(n):subid(n)+M;
            d = var(sid)';
            fd(:,n)=fft(d.*sw);
        end
        psd = sum(abs(fd(2:M/2+1,:)).^2,2)/N; %sum all ffts
        psd = psd/fs; %scaled power spectra
        f = linspace(0,fs/2,(M/2));
        alfa = 1 - 0.95;
        c = chi2inv([1-alfa/2 alfa/2],nu);
        c = nu./c; %spectra confidence intervals
        
        %save to structure
        if j==1,out.(fn{i}).f = f;end
        out.(fn{i}).psd(j,:) = psd;
        out.(fn{i}).ci(j,:) = c;
    end
end

%compute power spectra for h-bd and x-bd
for i = 1:2
    if i == 1
        var1 = data.h;
        
    else
        var1 = data.x;
    end
    var2 = data.bd;
    
    if length(var1) ~= length(var2)
        disp('h or x and bd must be the same length')
        return;
    end
    
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
        N = length(var);
        M = w;
        subid = [1 M/2:M/2:N]; %50% overlap
        sl = length(subid)-2;
        sw = hamming(M+1); %hamming window
        fd1 = zeros(M+1,sl); %preallocate
        fd2 = zeros(M+1,sl); %preallocate
        for n=1:sl
            sid = subid(n):subid(n)+M;
            d1 = v1(sid)';
            fd1(:,n)=fft(d1.*sw);
            d2 = v2(sid)';
            fd2(:,n)=fft(d2.*sw);
        end
        psd1 = fd1.*conj(fd1)/N^2;
        
        coh = abs(csd2).^2 ./(psd1.*psd2);
        phase = atan2(-real(csd),imag(csd))*180/pi;
    end
    
    
    



end