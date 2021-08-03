%Calculate Reynolds stress by calculating the cospectra of velocities, then
%filtering by the turbulent frequency band previously defined by our other
%analyses. We are adapting the CF method of Gerbi et al., 2008 here, but
%are not focusing on model fits to derive Reynolds stress. Instead, we
%simply calculate cospectra, define the cutoff frequency, and integrate the
%cospectra above this frequency.

clear
datdir = 'c:\Users\bkn5\Projects\Mekong_W2015\DataAnalysis\Paper1\';
fname = 'VTAvelocities.mat';
load([datdir fname])
fn = fieldnames(VTA);
ploton = 0;
cutoff = 6;
heading = 102;
%subdivide timeseries into sections 10 minutes long
intv = 10;
fs = 50;                                                            %Hz, this is the cutoff freq from the spectral slopes
avt = intv*fs*60;                                                   %samples/interval
for j = 2:4                                                         %loop by instrument
    x = VTA.(fn{j}).x;y = VTA.(fn{j}).y;z1 = VTA.(fn{j}).z1;z2 = VTA.(fn{j}).z2;
    
    %rotate to the cross-shore direction
    rot = heading*pi/180;
    x = x.*(ones(size(x))*cos(rot)) + ...
        y.*(ones(size(y))*sin(rot));
    y = -x.*(ones(size(x))*sin(rot)) + ...
        y.*(ones(size(y))*cos(rot));

    ind = [1 avt:avt:length(x)];
    for i = 1:length(ind)-1                                             %loop in time
        for ii = 1:35                                                   %loop by depth bin
            X = x(ind(i):ind(i+1),ii);
            Y = y(ind(i):ind(i+1),ii);
            Z1 = z1(ind(i):ind(i+1),ii);
            Z2 = z2(ind(i):ind(i+1),ii);
            if min(X) == 0 && max(X) == 0 || min(Y) == 0 && max(Y) == 0 %skip bins with zero velocity
                rss.(fn{j}).uw(i,ii) = 0;
                rss.(fn{j}).vw(i,ii) = 0;
                rss.(fn{j}).uwstar(i,ii) = 0;
                rss.(fn{j}).vwstar(i,ii) = 0;
            else
                %force inputs to be even
                if rem(length(X),2) == 1
                    X = X(1:end-1,:);
                end
                if rem(length(Y),2) == 1
                    Y = Y(1:end-1,:);
                end
                if rem(length(Z1),2) == 1
                    Z1 = Z1(1:end-1,:);
                end
                if rem(length(Z2),2) == 1
                    Z2 = Z2(1:end-1,:);
                end
                %compute along-beam velocities
                u = detrend(X);
                v = detrend(Y);
                w1 = detrend(Z1);
                w2 = detrend(Z2);
                theta = 30;                                                  %deg; half angle of beams to the vertical
                
                u1 = -u*sind(theta)-w1*cosd(theta);
                u2 = u*sind(theta)-w2*cosd(theta);
                u3 = -v*sind(theta)-w1*cosd(theta);
                u4 = v*sind(theta)-w2*cosd(theta);
                
                %compute velocity spectra
                [m,~] = size(u);
                nfft = 0.25*m;
                window = hanning(m,'periodic');
                noverlap = 0.5*length(window);
                
                [Su1,fu] = pwelch(u1,window,noverlap,nfft,fs,'onesided');
                [Su2,~] = pwelch(u2,window,noverlap,nfft,fs,'onesided');
                [Su3,~] = pwelch(u3,window,noverlap,nfft,fs,'onesided');
                [Su4,~] = pwelch(u4,window,noverlap,nfft,fs,'onesided');
                
                %calculate cospectra
                COuw = (Su1-Su2)./(4*cosd(theta)*sind(theta));
                COvw = (Su3-Su4)./(4*cosd(theta)*sind(theta));
                
                %by definition, the covariance of two signals is the integral of
                %the cospectrum
                uwstar = trapz(COuw);vwstar = trapz(COvw);
                
                %convert frequencies to wavenumbers via the Taylor Hypothesis
                kcu = cutoff./abs(nanmean(X));kcv = cutoff./abs(nanmean(Y)); %cutoff k
                ku = fu./abs(nanmean(X));kv = fu./abs(nanmean(Y));           %wavenumber (k)
                
                %compute Ogive curve (CDF) of cospectra
                oguw = cumtrapz(COuw);
                ogvw = cumtrapz(COvw);
                
                %turb energy is contained in low f in CDF functions
                uw = oguw(ku < kcu);
                vw = ogvw(kv < kcv);
                
                %peaks of the Ogive curves are the covariance
                %save variables in structure
                rss.(fn{j}).uw(i,ii) = mean(findpeaks(uw));
                rss.(fn{j}).vw(i,ii) = mean(findpeaks(vw));
                rss.(fn{j}).uwstar(i,ii) = uwstar;
                rss.(fn{j}).vwstar(i,ii) = vwstar;
            end
        end
    end
end
sfname = regexprep(fname,'velocities.mat','RS.mat');
save([datdir sfname],'rss')
if ploton
    for j = 2:4
        c = jet(35);
        figure
        for i = 1:35
            if any(rss.(fn{j}).uw(:,i) ~= 0)
                pf = polyfit(rss.(fn{j}).uw(:,i),rss.(fn{j}).uwstar(:,i),1);
                pv = polyval(pf,rss.(fn{j}).uw(:,i));
                resid = rss.(fn{j}).uw(:,i)-pv;
                SSresid = sum(resid.^2);
                SStotal = (length(rss.(fn{j}).uw(:,i))-1)*var(rss.(fn{j}).uw(:,i));
                rsq = 1-SSresid/SStotal;
                disp(['Bin # ' num2str(i) 'slope: ' num2str(pf(1)) ', r^2:' num2str(rsq)])
                plot(rss.(fn{j}).uw(:,i),rss.(fn{j}).uwstar(:,i),'o','color',c(i,:)), hold on
            end
        end
        pf = polyfit(rss.(fn{j}).uw(:,15),rss.(fn{j}).uwstar(:,15),1);
        pv = polyval(pf,rss.(fn{j}).uw(:,15));
        plot(rss.(fn{j}).uw(:,15),pv,'-k','LineWidth',1.5)
        xlabel('\int Co_u_w')
        ylabel('uw')
        
        figure
        for i = 1:35
            if any(rss.(fn{j}).vw(:,i) ~= 0)
                pf = polyfit(rss.(fn{j}).vw(:,i),rss.(fn{j}).vwstar(:,i),1);
                pv = polyval(pf,rss.(fn{j}).vw(:,i));
                resid = rss.(fn{j}).vw(:,i)-pv;
                SSresid = sum(resid.^2);
                SStotal = (length(rss.(fn{j}).vw(:,i))-1)*var(rss.(fn{j}).vw(:,i));
                rsq = 1-SSresid/SStotal;
                disp(['Bin # ' num2str(i) ' slope: ' num2str(pf(1)) ', r^2:' num2str(rsq)])
                plot(rss.(fn{j}).uw(:,i),rss.(fn{j}).uwstar(:,i),'o','color',c(i,:)), hold on
            end
        end
        pf = polyfit(rss.(fn{j}).vw(:,15),rss.(fn{j}).vwstar(:,15),1);
        pv = polyval(pf,rss.(fn{j}).vw(:,15));
        plot(rss.(fn{j}).vw(:,15),pv,'-k','LineWidth',1.5)
        xlabel('\int Co_v_w')
        ylabel('vw')
    end
end
