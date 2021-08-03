%simple Structure function method (to avoid having to use runVecPro.m)
%this code first runs through the instruments, then by time intervals
%defined in settings, then by the beam, then by the depth. TKE estimates
%are grown per beam per time interval until completion.
%%%%
clear
load('C:\Users\bkn5\Projects\Mekong_W2015\DataAnalysis\Paper1\HTAday1Vels.mat')
savedatdir = 'c:\Users\bkn5\Projects\Mekong_W2015\DataAnalysis\Turbulence\Paper1\';

fn = fieldnames(HTA);

%%%%
%Struct Fun basic settings:
window = 5; %5 minute averaging interval
step = 10; %10 minute step
fs = 50;
avt = step*fs*60; %samples/step
nwin = window*fs*60; %samples/window
cellsize = 0.9693;
r = (cellsize/1000)/cosd(30);                                               %converts vertical beam distance to along beam distance.
lags = 5;
TKE = struct();                                                             %initialize structures for data storage
Stat = struct();
%%%%
for i = 2:4                                                                 %loop through instruments VP1-VP3
    disp(['Analysing ' fn{i}])
    nsamp = length(HTA.(fn{i}).time);
    ind = [1 avt:avt:nsamp];
    [~,m] = size(HTA.(fn{i}).beam1);
    maxbin = m-lags;                                                        %TKE estimates are maximally length(maxbin) due to lag #
    %%%%    
    for ii = 1:length(ind)                                                  %loop through time, windowed to the settings
        if abs(nsamp-ind(ii)) < nwin                                        %skip the last few indexes approaching the end of the t-s
            continue
        else
            idx = ind(ii):ind(ii)+nwin-1;
            dat.beam1 = HTA.(fn{i}).beam1(idx,:);
            dat.beam2 = HTA.(fn{i}).beam2(idx,:);
            dat.beam3 = HTA.(fn{i}).beam3(idx,:);
            dat.beam4 = HTA.(fn{i}).beam4(idx,:);
            dfn = fieldnames(dat);
        end
        %%%%
        if ii == 1
            count = 1;
        elseif ii > 1
            count = count+lags;
        end
        for j = 1:4                                                         %loop through beams
            %Define variables:
            itt = 1;                                                        %iteration number
            E = zeros(maxbin,1);
            D = zeros(maxbin,lags);
            pf = zeros(maxbin,2);
            pv = zeros(maxbin,lags);
            pf2 = zeros(maxbin,2);
            pv2 = zeros(maxbin,lags);
            pval = zeros(maxbin,1);
            rsq = zeros(maxbin,1);
            rsq2 = zeros(maxbin,1);
            bfit = zeros(maxbin,1);
            R = zeros(maxbin,lags);
            rr = zeros(1,lags);
            N = zeros(1,maxbin);
            epsilon = zeros(1,maxbin);
            epserr = zeros(1,maxbin);
            %%%%
            beam = detrend(dat.(dfn{j}),'constant');                        
            while itt <= maxbin                                             %loop in depth
                d = zeros(length(idx),lags);
                idl = itt:itt+lags;
                idt = 1:length(idx);
                for jj = 1:lags                                             %compute velocity differences (vel(z,1) - vel(z,2:lags))^2
                    d(:,jj) = (beam(idt,itt)-(beam(idt,idl(jj+1)))).^2;
                    D(itt,jj) = nanmean(d(:,jj));
                    rr(:,jj) = r*jj;
                end
                R(itt,:) = rr.^(2/3);                                       %along beam distance r^2/3
                [pf(itt,:),S] = polyfit(R(itt,:),D(itt,:),1);               %linear regression of D along r^2/3
                [~,~,~,~,stats] = regress(D(itt,:)',[ones(size(R(itt,:)')) R(itt,:)'],0.05);
                pval(itt,:) = stats(3);                                     %p-value with 95% confidence intervals
                pv(itt,:) = polyval(pf(itt,:),R(itt,:));
                CI = polyparci(pf(itt,:),S,0.95);                           %95% confidence interval on r^2/3
                CImax = max(CI(:,1));
                A = pf(itt,1);
                N(itt) = pf(itt,2);
                epsilon(itt) = (A/2.1)^(3/2);                               %units of W/m^3
                epsmax = (CImax/2.1)^(3/2);
                epserr(itt) = epsmax-epsilon(itt);
                %%%%                                                        %calculate residuals as a measure of fit
                resid = D(itt,:)-pv(itt,:);                                 
                SSresid = sum(resid.^2);
                SStotal = (length(D(itt,:))-1)*var(D(itt,:));
                rsq(itt,:) = 1-SSresid/SStotal;
                %%%%                                                        %iteratively calculate a best-fit slope from range 0.1-1.5 (2/3 is 0.6667)
                b = 0.1:0.01:1.5;
                rsqb = zeros(length(b),1);
                pfb = zeros(length(b),2);
                pvb = zeros(length(b),lags);
                for jj = 1:length(b)
                    Rb = rr.^b(jj);
                    pfb(jj,:) = polyfit(Rb,D(itt,:),1);                      
                    pvb(jj,:) = polyval(pfb(jj,:),Rb);
                    resid = D(itt,:)-pvb(jj,:);                                   
                    SSresid = sum(resid.^2);
                    SStotal = (length(D(itt,:))-1)*var(D(itt,:));
                    rsqb(jj) = 1-SSresid/SStotal;
                end
                [rsq2(itt,:),id] = max(rsqb);
                pf2(itt,:) = pfb(id,:);
                pv2(itt,:) = pvb(id,:);
                bfit(itt,:) = b(id);
                %%%%
                itt = itt+1;
            end
            E(1:length(epsilon),1) = epsilon';                              %TKE dissipation rate
            E(abs(imag(E))>0) = NaN;                                        %filter non-real numbers from E estimates, Noise
            E(abs(imag(N'))>0) = NaN;
            epserr(abs(imag(epserr))>0) = NaN;                              %do the same for error estimates on E
            rsq(isinf(rsq)) = NaN;
            rsq2(isinf(rsq)) = NaN;
            
            %%%%
            %Save variables:
            TKE.(dfn{j}).E(1:maxbin,ii) = E;                                %save epsilon out to structure TKE by beam
            TKE.(dfn{j}).Emaxerr(1:maxbin,ii) = epserr;                     %maximum error in epsilon estimates from r^2/3 fits
            TKE.(dfn{j}).R(1:lags,ii) = R(1,:);                             %R is the along-beam distances from the r^2/3 fits
            TKE.(dfn{j}).r(1:lags,ii) = rr';                                %r is the along-beam distances from the r fits
            TKE.(dfn{j}).A(1:maxbin,ii) = pf(:,1);                          %A is the slope from the r^2/3 fits
            TKE.(dfn{j}).a(1:maxbin,ii) = pf2(:,1);                         %a is the slope from the ln(r) fits
            TKE.(dfn{j}).bfit(1:maxbin,ii) = bfit(:,1);                     %bfit is the slope [0.1-1.5] that best fits the data
            TKE.(dfn{j}).RSQ(1:maxbin,ii) = rsq;                            %RSQ is r-squared from the r^2/3 fits
            TKE.(dfn{j}).rsq(1:maxbin,ii) = rsq2;                           %rsq is r-squared from the ln(r) fits
            TKE.(dfn{j}).pval(1:maxbin,ii) = pval;                          %p-value statistics from 2/3 fits
            if ii == 1;
                TKE.(dfn{j}).Yfit(:,1:lags) = pv;
                TKE.(dfn{j}).yfit(:,1:lags) = pv2;
                TKE.(dfn{j}).D(:,1:lags) = D;                               %D has dimensions: bins x lags. When saved into the structure, D is bins x lags*length(ind)
            elseif ii > 1;
                TKE.(dfn{j}).Yfit(:,count:count+lags-1) = pv;
                TKE.(dfn{j}).yfit(:,count:count+lags-1) = pv2;
                TKE.(dfn{j}).D(:,count:count+lags-1) = D;                   %Each step in time is D(n:n+5) where n is the timestep. 
            end
            TKE.(dfn{j}).lags = lags;
        end
    end
    Stat.(fn{i}) = TKE;                                                     %save variables out to main structure
end
disp('Finished!')

disp('Saving HTA 07/03/2015 TKE file')
save([savedatdir 'HTAday1TKE'],'Stat','-v7.3')

