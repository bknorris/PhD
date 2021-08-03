%Plot Ub & BLE (bed movement) for an experiment.
%
% Paper 3: Sediment Motion in Mangroves
%
% This is Version 1.0 of this script.
%
% This script tries to correlate bed movement from the VP with near-bottom
% orbital velocities calculated from VP velocities. Also, PSD and CPSD are
% attempted with the bed (again).
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
bdfiles = {'VTA1a_bdtrace.mat'};
dir2 = 'd:\Projects\Mekong_W2015\DataAnalysis\Paper3\VPs\';
vfiles = {'13March2015a_Vels.mat'};
%
%Pressure sensor
load('D:\Projects\Mekong_W2015\Data\RBR\Duet\DPS2\Duet_140315.mat')
ll = sin(RBR.Metadata.latitude/57.29578)^2;
zp = str2double(regexp(RBR.Metadata.hab,'^\S+','match'))/1000;
%
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
t1 = {'13-Mar-2015 16:41:31'};
t2 = {'13-Mar-2015 17:07:55'};
for i = 1
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
        %Interpolate x,y (downsample) to 10Hz
        x10 = zeros(length(bdt),35);
        y10 = zeros(length(bdt),35);
        z10 = zeros(length(bdt),35);
        %Find unique values to avoid error with interp1
        bu = unique(bdt);
        [vu,id] = unique(vpt);
        for j = 1:35
            x10(:,j) = interp1(vu,x(id,j),bu);
            y10(:,j) = interp1(vu,y(id,j),bu);
            z10(:,j) = interp1(vu,y(id,j),bu);
        end
        %Prepare pressure data
        cp3 = find(RBR.Datetime >= start & RBR.Datetime <= stop);
        pt = RBR.Datetime(cp3);
        pres = RBR.SeaPres(cp3);
        p10 = interp1(pt,pres,bu);
        %Compute spectra:
        fs = 10;
        avt = 10*fs; %seconds
        nwin = 120*fs; %seconds
        lt = length(x10);
        ind = [1 avt:avt:lt]; 
        %Define variables
        ub = zeros(floor(lt/avt),1);
        up = zeros(floor(lt/avt),1);
        wp = zeros(floor(lt/avt),1);
        uw = zeros(floor(lt/avt),1);
        bm = zeros(floor(lt/avt),1);
        bvar = zeros(floor(lt/avt),1); 
        xbh = zeros(floor(lt/avt),301);
        E = zeros(floor(lt/avt),1);
        t1 = zeros(floor(lt/avt),1);
        for j = 1:floor(lt/avt)
            if abs(lt-ind(j)) < nwin  %skip the last few indexes approaching the end of the t-s
                continue
            else
                idx = ind(j):ind(j)+nwin-1;
                t1(j) = bdt(idx(1));
                xm = detrend(nanmean(x10(idx,bins),2));
                ym = detrend(nanmean(y10(idx,bins),2));
                zm = detrend(nanmean(z10(idx,bins),2));
                bh = detrend(bdh(idx));
                P = p10(idx);pm = detrend(P);
                g = zeros(length(P),1);h = zeros(length(P),1);
                for jj = 1:length(P)
                    g(jj,:) = 9.780318*(1.0+(5.2788E-3+2.36E-5*ll)*ll)+1.092E-6*P(jj,:);
                    h(jj,:) = ((((-1.82E-15*P(jj,:)+2.279E-10)*P(jj,:)-2.2512E-5)*P(jj,:)+9.72659)*P(jj,:))/g(jj,:);
                end
                h = mean(h);
                g = mean(g);
                %Spectra
                M = floor(nwin/12);
                [Cxx,f] = pwelch(xm,hamming(M),[],nwin/2,fs);
                [Cyy,~] = pwelch(ym,hamming(M),[],nwin/2,fs);
                [Czz,~] = pwelch(zm,hamming(M),[],nwin/2,fs);
                [Cpp,~] = pwelch(pm,hamming(M),[],nwin/2,fs);  
                %compute wave-cutoff (zid)
                dof = 2*floor(nwin/M);ci = cohere_signif_level(dof);
                [msc,~] = mscohere(xm,zm,hamming(floor(M)),[],nwin/2,fs);
                wband = find(f>0.05 & f<1.2);
                zid = wband((msc(wband) > ci));
                Guv = Cxx(zid)+Cyy(zid);
                ub(j) = sqrt(2*sum(Guv)*df); %orbital velocity
                [xbh(j,:),~] = mscohere(xm,bh,hamming(floor(M)),[],nwin/2,fs);
                %Following Bricker & Monismith, 2007:
                %Integrate PSD of X (Y,Z) and subtract from integral of X
                %(Y,Z) for values of CSD(XZ) > ci. This is the total stress
                %below the wave band. The PSD of X (Y,Z) for CSD(XZ) > ci is
                %the wave band, from which orbital velocities may be calculated.
                Cxx(zid(1:end)) = NaN;Cxx = fixgaps(Cxx);
                Czz(zid(1:end)) = NaN;Czz = fixgaps(Czz);
                
                up(j) = trapz(Cxx); %"u prime"
                wp(j) = trapz(Czz); %"w prime"
                
                uw(j) = up(j)*wp(j); %Cross-shore, vertical Reynolds stress
                
                %Compute Hs, Wave energy from pressure t-s
                df = f(3)-f(2);
                omega = 2*pi.*f;
                k = qkhf(omega,h)./h;
                kh = k*h;
                kz = k*zp;
                attn = cosh(kz)./cosh(kh);
                attn(attn<0.2) = 0.2;
                Spp = Cpp./(attn.^2);           %surface elevation spectrum
                m0 = sum(Spp(zid)*df);
                Hs = 4*sqrt(m0);                %sig. wave height
                rho = 1010; %kg/m^3
                E(j) = (1/8)*rho.*g.*(Hs/sqrt(2)).^2; %wave energy
                %Bottom movement
                bm(j) = bdh(idx(1));
                bvar(j) = mean(gradient(bdh(idx)));
            end
        end

        up(bm==0)=[];wp(bm==0)=[];uw(bm==0)=[];t1(bm==0)=[];E(bm==0)=[];
        bvar(bm==0)=[];ub(bm==0)=[];bm(bm==0)=[];
        figure
        sp(1) = subplot(311);
        plot(t1,bm,'k','linewidth',1.5)
        title('BLE')
        ylabel('m')
        sp(2) = subplot(312);
        plot(t1,E,'k','linewidth',1.5)
        title('Wave Energy')
        ylabel('m^2')
        sp(3) = subplot(313);
        plot(t1,uw,'k','linewidth',1.5)
        title('u_b')
        ylabel('m/s')
        set([sp(1) sp(2)],'xticklabel',[])
        set(sp,'xlim',[t1(1) t1(end)])
        datetickzoom('x','HH:MM','keepticks','keeplimits')
        linkaxes(sp,'x')

        
        
        figure
        set(gcf,'Renderer','zbuffer')
        sbmin = min(bvar);sbmax = max(bvar);
        cax = linspace(sbmin,sbmax,length(bvar));
        cl = brewermap(length(cax),'RdBu');cl = flipud(cl);
        pf = polyfit(bm,uw,1);
        xs = linspace(min(bm),max(bm),length(bm));
        pv = polyval(pf,xs);
        for j = 1:length(bm)
            tmp = abs(cax-bvar(j));
            [~,I]=min(tmp);
            plot(bm(j),uw(j),'o',...
                'markersize',8,...
                'markerfacecolor',cl(I,:),...
                'markeredgecolor','k'),hold on
        end
        plot(xs,pv,'-k','linewidth',1.5)
        colormap(cl)
        cb = colorbar('v');
        caxis([sbmin sbmax])
        ylabel(cb,'$log_{10}\left (\sigma_{b,int}^{2}\right )$',...
            'interpreter','latex')
        ylabel('$\left | \epsilon \right | (Wkg^{-1})$','interpreter','latex')
        xlabel('$BLE (m)$','interpreter','latex')
        set(gca,'position',[0.1 0.13 0.65 0.8],...
            'xlim',[min(bm) max(bm)])
        prettyfigures('font','Cambria','text',12,'labels',14,...
            'box',1)
    end
end

            
            
