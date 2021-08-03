%Plot TKE & BLE (bed movement) for an experiment. Then, calculate
%xcorr(TKE,BLE) if sensible.
%
% Paper 3: Sediment Motion in Mangroves
%
% This is Version 1.2 of this script
%
% V 1.0: xcorr of TKE and BLE was inconclusive. Sometimes TKE preceeded
% BLE, othertimes the converse.
% V 1.1: Tried calculating spectra with FFT. Required collocated pressure
% sensor to find fwc (wave-band cutoff f). New method (from Bricker &
% Montismith, 2007) uses CSD of u and w to find fwc. 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear
close all
expname = 'VTA';
start = datenum('13-Mar-2015 00:00:00');
stop = datenum('15-Mar-2015 00:00:00');
wbbl = [0.0025 0.0025]; %from WaveBoundaryThickness.m
bins = 13:17;
hoffset = 99; %difference between vp heading and transect (deg)
%
dir1 = 'd:\Projects\Mekong_W2015\DataAnalysis\Paper3\BottomTrack\';
bdfiles = {'VTA1a_bdtrace.mat';'VTA2b_bdtrace.mat'};
dir2 = 'd:\Projects\Mekong_W2015\DataAnalysis\Paper3\VPs\';
vfiles = {'13March2015a_Vels.mat';'14March2015b_Vels.mat'};
%
%Pressure sensor
load('D:\Projects\Mekong_W2015\Data\RBR\Duet\DPS2\Duet_140315.mat')
ll = sin(RBR.Metadata.latitude/57.29578)^2;
zp = str2double(regexp(RBR.Metadata.hab,'^\S+','match'))/1000;
%Plot colors
if strfind(expname,'VTA')
    cl = [0 0 0];
else
    cl = brewermap(5,'Blues');
    cl = cl(3:5,:);
end
%Rangebins height
rb = [0.0400,0.0410,0.0419,0.0429,0.0439,...
0.0448,0.0458,0.0467,0.0477,0.0487,...
0.0496,0.0506,0.0516,0.0525,0.0535,...
0.0545,0.0554,0.0564,0.0574,0.0583,...
0.0593,0.0602,0.0612,0.0622,0.0631,...
0.0641,0.0651,0.0660,0.0670,0.0680,...
0.0689,0.0699,0.0708,0.0718,0.0728];
sdir = 'd:\Projects\Mekong_W2015\Figures\Paper3\BottomTracking\';

%Hone in on a specific time period when TKE is high; look for xcorrelation
%between the two signals
t1 = {'13-Mar-2015 16:41:31';'14-Mar-2015 13:06:52'};
t2 = {'13-Mar-2015 17:07:55';'14-Mar-2015 13:41:48'};
for i = 1:length(bdfiles)
    start = datenum(t1{i});
    stop = datenum(t2{i});
    %Load the data
    bd = load([dir1 bdfiles{i}]);
    fn = whos('-file',[dir2 vfiles{i}]);fn = {fn.name};
    v = matfile([dir2 vfiles{i}]);
    if strfind(expname,'VTA')
        ni = 1;
    else
        ni = 1:3;
    end
    for ii = ni
        %Prepare bottom tracking data
        cp1 = find(bd.(fn{ii}).time >= start & bd.(fn{ii}).time <= stop);
        bdt = bd.(fn{ii}).time(cp1);
        bdh = bd.(fn{ii}).bdist(cp1);
        bid = find(~isnan(bdh),1,'first');bdi = bdh(bid);
        bdh = bdi-bdh;
        vph = bdi-rb;
        %Prepare velocity data
        vp = v.(fn{ii});
        cp2 = find(vp.time >= start & vp.time <= stop);
        vpt = vp.time(cp2);
        x = vp.x(cp2,:);
        y = vp.y(cp2,:);
        z = (vp.z1(cp2,:)+vp.z2(cp2,:))./2;
        %Rotate horiz. velocities to along & cross shore
        rot = hoffset(ii)*pi/180;
        x = x.*(ones(size(x))*cos(rot)) + ...
            y.*(ones(size(y))*sin(rot));
        y = -y.*(ones(size(y))*sin(rot)) + ...
            x.*(ones(size(x))*cos(rot));
        %Prepare pressure data
        cp3 = find(RBR.Datetime >= start & RBR.Datetime <= stop);
        pt = RBR.Datetime(cp3);
        pres = RBR.SeaPres(cp3);
        p50 = interp1(pt,pres,vpt);
        %Compute turbulence statistics:
        %Spectral settings
        fs = 50;
        avt = 5; %to get to 10Hz
        nwin = 60*fs;
        lt = length(x);
        ind = [1 avt:avt:lt]; 
        %Define variables
        up = zeros(floor(lt/avt),1);
        vp = zeros(floor(lt/avt),1);
        wp = zeros(floor(lt/avt),1);
        uw = zeros(floor(lt/avt),1);
        vw = zeros(floor(lt/avt),1);
        k = zeros(floor(lt/avt),1);
        ub = zeros(floor(lt/avt),1);
        matlabpool(4)
        parfor j = 1:floor(lt/avt)
            if abs(lt-ind(j)) < nwin  %skip the last few indexes approaching the end of the t-s
                continue
            else
            idx = ind(j):ind(j)+nwin-1;
            end
            %
            X = detrend(nanmean(x(idx,bins),2));
            Y = detrend(nanmean(y(idx,bins),2));
            Z = detrend(nanmean(z(idx,bins),2));
            P = p50(idx);
            %Compute h - water depth
            g = zeros(length(P),1);h = zeros(length(P),1);
            for jj = 1:length(P)
                g(jj,:) = 9.780318*(1.0+(5.2788E-3+2.36E-5*ll)*ll)+1.092E-6*P(jj,:);
                h(jj,:) = ((((-1.82E-15*P(jj,:)+2.279E-10)*P(jj,:)-2.2512E-5)*P(jj,:)+9.72659)*P(jj,:))/g(jj,:);
            end
            h = mean(h);
            g = mean(g);
            P = detrend(P);
            M = max(size(X))/8;
            %Spectra settings:
            [Cxx,f] = pwelch(X,hamming(floor(M)),[],nwin/2,fs);
            df = f(3)-f(2);
            [Cyy,~] = pwelch(Y,hamming(floor(M)),[],nwin/2,fs);
            [Czz,~] = pwelch(Z,hamming(floor(M)),[],nwin/2,fs);
            wband = find(f>0.05 & f<1.2);
            dof = 2*(nwin/M);ci = cohere_signif_level(dof);
            [msc,~] = mscohere(X,Z,hamming(floor(M)),[],nwin/2,fs);
            xzid = wband((msc(wband) > ci));
            [msc,~] = mscohere(Y,Z,hamming(floor(M)),[],nwin/2,fs);
            yzid = wband((msc(wband) > ci));
            
            %Following Bricker & Monismith, 2007: 
            %Integrate PSD of X (Y,Z) and subtract from integral of X
            %(Y,Z) for values of CSD(XZ) > ci. This is the total stress 
            %below the wave band. The PSD of X (Y,Z) for CSD(XZ) > ci is
            %the wave band, from which orbital velocities may be calculated.
            Cxx(xzid(2:end)) = NaN;Cxx = fixgaps(Cxx);
            Cyy(yzid(2:end)) = NaN;Cyy = fixgaps(Cyy);
            Czz(xzid(2:end)) = NaN;Czz = fixgaps(Czz);
            
            up(j) = trapz(Cxx); %"u prime"
            vp(j) = trapz(Cyy); %"v prime"
            wp(j) = trapz(Czz); %"w prime"
            
            uw(j) = up(j)*wp(j); %Cross-shore, vertical Reynolds stress
            vw(j) = vp(j)*wp(j); %Along-shore, " "
            
            k(j) = 0.5*(up(j)^2+vp(j)^2+wp(j)^2); %TKE
            if length(xzid) > length(yzid)
                zid = yzid;
            elseif length(xzid) < length(yzid)
                zid = xzid;
            else
                zid = xzid;
            end
            Guv = Cxx(zid)+Cyy(zid);
            ub(j) = sqrt(2*sum(Guv)*df);   
        end
        matlabpool CLOSE
        %Compure xcorr
        figure
        fs = 10;
        avt = fs*60; %60 second
        nwin = fs*120;
        nsamp = length(ub);
        lt = length(bdh);
        ind = [1 avt:avt:lt];
        cl = brewermap(floor(lt/avt),'Blues');
        tD = zeros(floor(lt/avt),1);
        tt = zeros(floor(lt/avt),1);
        for j = 1:length(ind)
            if abs(nsamp-ind(j)) < nwin
                continue
            else
                idx = (j-1)*avt+1:j*avt;
                [acor,lag]=xcorr(ub(idx),bdh(idx));
                [~,I] = max(abs(acor));
                tD(j) = lag(I)/fs;
                tt(j) = bdt(idx(1));
                plot(lag/fs,acor,'color',cl(j,:)),hold on
            end
        end
        tD(tt==0)=[];tt(tt==0)=[];
        xlabel('Lag time (s)')
        ylabel('Correlations R_u_b_,_B_L_E')
        figure
        plot(tt,tD,'k')
        set(gca,'xlim',[tt(1) tt(end)])
        datetickzoom('x','HH:MM','keepticks','keeplimits')
        xlabel(['Time on ' datestr(tt(1),'dd-mm-yy')])
        ylabel('Time lag (s)')
        title('Time lag from xcorr(ub,BLE)')
%         export_fig(figure(2),[sdir expname '_ubBLE_timelag'],'-png')
%         export_fig(figure(6),[sdir expname '_ubBLE'],'-png')
%         PlotBLExcorr
    end
end
%             %Interpolate to 10Hz t-s (downsample vels & pres)
%         p10 = interp1(pt,pres,bdt);
%         x10 = zeros(length(bdt),min(size(x)));
%         y10 = zeros(length(bdt),min(size(y)));
%         z10 = zeros(length(bdt),min(size(z)));
%         for k = 1:min(size(x))
%             x10(:,k) = interp1(vpt,x(:,k),bdt);
%             y10(:,k) = interp1(vpt,y(:,k),bdt);
%             z10(:,k) = interp1(vpt,z(:,k),bdt);
%         end
%         %Compure xcorr
%         fs = 10;
%         avt = fs*60; %60 second
%         lt = length(p10);
%         ind = [1 avt:avt:lt];
%         for k = 1:floor(lt/avt)
%             %Average x,y,z across bins
%             X = nanmean(x10((k-1)*avt+1:k*avt,bins),2);
%             Y = nanmean(y10((k-1)*avt+1:k*avt,bins),2);
%             Z = nanmean(z10((k-1)*avt+1:k*avt,bins),2);
%             P = p10((k-1)*avt+1:k*avt);
%             BLE = bdh((k-1)*avt+1:k*avt);
%             U = sqrt(X.^2+Y.^2);
%             W = Z;
%             %Calculate u' and w'
%             X = detrend(X);
%             Y = detrend(Y);
%             Z = detrend(Z);
%             P = detrend(P);
%             %Spectra settings:
%             M = 60;
%             sw = hamming(M);
%             [Cxx,f]=pwelch(X);
            
            
            
