clear
%This script identifies strongly onshore and offshore flow near the bed,
%then calculates the dissipation rate from these samples, and plots a
%histogram of the occurences for the paper.
%

%load vectrino/turbulence files from HTA1
datdir = 'g:\Mekong_W2015\DataAnalysis\Paper3\VPs\';
m = matfile([datdir '7March2015_Vels.mat']);
fn = whos(m);fn = {fn.name};
start = datenum(2015,03,07,14,00,00);
stop = datenum(2015,03,07,17,00,00);
%Find onshore/offshore flow
bin = 5;
fs = 50;
data = struct();
for i = 1:3
    dat = m.(fn{i});
    t = dat.time;
    id = find(t>=start & t<= stop);
    x = mean(dat.x(id,1:bin),2);y = mean(dat.y(id,1:bin),2);
    rot = (pi*340)/180;
    T = [cos(rot) -sin(rot);...
        sin(rot) cos(rot)];
    vels = [x y];
    V = vels*T';
    x = V(:,1);y = V(:,2);
    onshore = find(x>0.3);
    offshore = find(x<-0.3);
    %save to structure for plotting
    data.(fn{i}).x = x;
    data.(fn{i}).onshore = onshore;
    data.(fn{i}).offshore = offshore;
    %calculate structure function dissipation rate
    cellsize = 0.9693;
    r = (cellsize/1000)/cosd(30);   %converts vertical beam distance to along beam distance.                                                                    %vertical velocities are not converted.
    lags = 5;
    ids = {offshore; onshore};
    maxbin = 35-lags;
    eps = zeros(2,2);
    for k = 1:2
        idx = ids{k};
        DATA.z1 = dat.z1(idx,:);  %we are now using vertical beams only so Cv2 = 2.1.
        DATA.z2 = dat.z2(idx,:);
        dfn = fieldnames(DATA); 
        for j = 1:length(dfn)   %loop through beams
            %Define variables:
            itt = 1;             %iteration number
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
            beam = detrend(DATA.(dfn{j}),'constant');
            while itt <= maxbin                        %loop in depth
                d = zeros(length(idx),lags);
                idl = itt:itt+lags;
                idt = 1:length(idx);
                for jj = 1:lags                         %compute velocity differences (vel(z,1) - vel(z,2:lags))^2
                    d(:,jj) = (beam(idt,itt)-(beam(idt,idl(jj+1)))).^2;
                    D(itt,jj) = nanmean(d(:,jj));
                    rr(:,jj) = r*jj;
                end
                R(itt,:) = rr.^(2/3);                    %along beam distance r^2/3
                [pf(itt,:),S] = polyfit(R(itt,:),D(itt,:),1);  %linear regression of D along r^2/3
                A = pf(itt,1);
                N(itt) = pf(itt,2);
                epsilon(itt) = (A/2)^(3/2);              %units of W/m^3
                itt = itt+1;
            end
            E(1:length(epsilon),1) = epsilon';            %TKE dissipation rate
            E(abs(imag(E))>0) = NaN;                      %filter non-real numbers from E estimates, Noise
            E(abs(imag(N'))>0) = NaN;
            eps(k,j) = nanmean(E(1:bin));
        end
    end
    eps = mean(eps);
    data.(fn{i}).ton = eps(1);
    data.(fn{i}).toff = eps(2);
end
%Plot Routine
f1 = figure(1);
set(f1,'PaperOrientation','portrait',...
    'position',[400 200   500   400]);
set(gcf,'color','w','paperpositionmode','auto')
ct = [1 4 7];
vp = [1 2 3];
c = flipud([0.1 0.1 0.1;0.5 0.5 0.5;0.7 0.7 0.7]);
for i = 1:3
    cts = ct(i);
    off = length(data.(fn{vp(i)}).offshore);
    on = length(data.(fn{vp(i)}).onshore);
    center = [cts-0.5 cts+0.5];
    br = bar(center,[off on]); hold on
    set(br,'facecolor',c(i,:),...
        'linewidth',1.5)
%     text(cts-0.7,240,'Off     On')
end