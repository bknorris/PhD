%Calculate power spectra of depth (ADCP), velocity (VP) and bed trace (VP).
%Additionally, compute the CPSD of velocity and bed trace (see Puleo et al.
%2014).
%
% Paper 3: Sediment Motion in Mangroves
%
% This is Version 1.0 of this script
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

F2F2 = struct();
ind = [1 avt:avt:nsamp];
for j = 1:length(ind)
    if abs(nsamp-ind(j)) < nwin  %skip the last few indexes approaching the end of the t-s
        continue
    else
        idx = ind(j):ind(j)+nwin-1;
    end
    X = x(idx);
    Y = y(idx);
    H = depth(idx);
    BLE = bd25(idx);

    %detrend
    X = detrend(X);Y = detrend(Y);
    H = detrend(H);BLE = detrend(BLE);
    
    w = hanning(floor(max(size(H))/16));
    [Chh,F] = pwelch(H,w,0,[],fs);
    nu = 2*floor(length(H)/length(w));
    err_low = nu/chi2inv(.05/2,nu);
    err_high = nu/chi2inv(1-.05/2,nu);

    [Cxx,~] = pwelch(X,hanning(swin),[],nwin/5,fs,'power','ConfidenceLevel',0.95);
    [Cyy,~] = pwelch(Y,hanning(swin),[],nwin/5,fs,'power','ConfidenceLevel',0.95);
    [Cble,~] = pwelch(BLE,hanning(swin),[],nwin/5,fs,'power','ConfidenceLevel',0.95);
    
    F2F2.(fn1{ii}).time(j) = tb(idx(1));
    F2F2.(fn1{ii}).F = F;
    F2F2.(fn1{ii}).Chh(j,:) = Chh;
    F2F2.(fn1{ii}).Cxx(j,:) = Cxx;
    F2F2.(fn1{ii}).Cyy(j,:) = Cyy;
    F2F2.(fn1{ii}).Cble(j,:) = Cble;
    
    %cross spectra, MSC and phase
    [Cxb,~] = cpsd(X,BLE,hanning(swin),[],nwin,fs);
    [MSC,~] = mscohere(X,BLE,hanning(swin),[],nwin,fs);
    phase = rad2deg(angle(Cxb));
    
    F2F2.(fn1{ii}).Cxb(j,:) = Cxb;
    F2F2.(fn1{ii}).MSCxb(j,:) = MSC;
    F2F2.(fn1{ii}).Pxb(j,:) = phase;
    
    [Chb,F] = cpsd(H,BLE,hanning(swin),[],nwin,fs);
    [MSC,~] = mscohere(H,BLE,hanning(swin),[],nwin,fs);
    phase = rad2deg(angle(Chb));
    
    F2F2.(fn1{ii}).Chb(j,:) = Chb;
    F2F2.(fn1{ii}).MSChb(j,:) = MSC;
    F2F2.(fn1{ii}).Phb(j,:) = phase;
end
