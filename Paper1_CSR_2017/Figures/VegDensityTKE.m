%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%We're going to plot TKE dissipation rates (averaged over 10 minute
%windows) against the vegetation statistics (phi) to see if there is a
%trend in TKE with increasing vegetation density.
clear
%Deployment orientations (2014):
%Q2 - MA1; Vertical VPs, VP1 east and VP3 west; HAB 139/141/138 profiles between 99 and 69mm
%Q3a - FSS.1; Vertical VPs, VP1 shoreward VP3 landward; HAB 81/83/77 profiles between 41 and 11mm
%Q3b - FSS.2; Horiz VPs, VP1 bottom VP3 top; HAB 83/179/279; profiles at 83mm (minimum)
%Q0 - Control (no quadrat); Horiz VPs, VP1 bottom VP3 top; HAB 101/204/300; profiles at 101mm (minimum)
%Q4 - FSS2; Horiz VPs, VP1 bottom VP3 top; HAB 84/174/272; profiles at 84mm (minimum)
%Q5 - F2F (fringe, VP2), Vertical VPs; HAB 150/82/150; profiles between 42 and 12mm
%Q6 - F2F (forest, VP3), Vertical VPs; HAB 150/82/150; profiles between 110 and 80mm
%Q7 - DPS, Horiz VPs, VP1 bottom VP3 top; HAB 108/208/308; profiles at 108mm (minimum)
%Deployment orientations (2015):
%Q1_1A & Q1_4A F2F2; VPs 2 and 3, both days. VP1 has no vegetation.
%Q2 - FSS3.1;
%Q2 - FSS3.2;
%Q2 - FSS3.3;
%Q3A & Q3B  - F2F3; VPs 2 and 3, both days. VP1 has no vegetation.
%Q4 - DPS2;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%BEGIN DATA CONGLOMERATOR%%%%
z = [99 101 98;41 43 37;83 179 279;101 204 300;84 174 272;NaN 42 NaN;NaN NaN 110;NaN 42 NaN;NaN NaN 110;108 208 308];

%%%%Load 2014 Data
veldir = 'd:\Projects\Mekong_F2014\DataAnalysis\Paper2\AveragedVelocities\';
wvdir = 'd:\Projects\Mekong_F2014\DataAnalysis\Paper2\Spectra\';
tkedir = 'd:\Projects\Mekong_F2014\DataAnalysis\Paper2\TKE\';
vegdir = 'd:\Projects\Mekong_F2014\DataAnalysis\Paper2\VegStats\';
velfiles = dir(veldir);velfiles = {velfiles.name};
tkefiles = dir(tkedir);tkefiles = {tkefiles.name};
vegfiles = dir(vegdir);vegfiles = {vegfiles.name};
wvfiles = dir(wvdir);wvfiles = {wvfiles.name};

%extract TKE estimates from TKE files
horiz = [1 1 0 0 1 0 1 0];
eps = zeros(length(tkefiles)-2,3);
wvrms = zeros(length(tkefiles)-2,1);
stdev = zeros(length(tkefiles)-2,3);
c = 1; %index counter
for i = 3:length(tkefiles)
    load([tkedir tkefiles{i}])
    load([veldir velfiles{i}])
    load([wvdir wvfiles{i}])
    %crop wavestats to the TKE deployment record: from closest values
    vl = length(Stat.vpro1.time);
            ind = zeros(1,vl);
            for j = 1:vl 
                tmp = abs(wvstats.time-Stat.vpro1.time(j));
                [~,idx] = min(tmp);
                ind(j) = idx;
            end
        wvel = wvstats.rorb(ind);
        wrms = sqrt(sum(wvel.^2)/numel(wvel));
    if horiz(c) == 1  %if the instrument is horizontal
        tk = (Stat.vpro1.beam1.E(15,1:vl)+Stat.vpro1.beam3.E(15,1:vl))./2;
        av = Avgs.vpro1.x(1,1:vl)+wrms;                                     %normalize by mean(u)+wrms (cubed); gives units of dissipation rate per meter
        ntke = tk./(av.^3);stdv(1) = std(ntke);
        e(1) = nanmean(ntke);
        tk = (Stat.vpro2.beam1.E(15,1:vl)+Stat.vpro2.beam3.E(15,1:vl))./2;
        av = Avgs.vpro2.x(1,1:vl)+wrms;
        ntke = tk./(av.^3);stdv(2) = std(ntke);
        e(2) = nanmean(ntke);
        tk = (Stat.vpro3.beam1.E(15,1:vl)+Stat.vpro3.beam3.E(15,1:vl))./2;
        av = Avgs.vpro3.x(1,1:vl)+wrms;
        ntke = tk./(av.^3);stdv(3) = std(ntke);
        e(3) = nanmean(ntke);
    elseif horiz(c) == 0
        tk = (Stat.vpro1.beam1.E(1,1:vl)+Stat.vpro1.beam3.E(1,1:vl))./2;
        av = Avgs.vpro1.x(1,1:vl)+wrms;
        ntke = tk./(av.^3);stdv(1) = std(ntke);
        e(1) = nanmean(ntke);
        tk = (Stat.vpro2.beam1.E(1,1:vl)+Stat.vpro2.beam3.E(1,1:vl))./2;
        av = Avgs.vpro1.x(1,1:vl)+wrms;
        ntke = tk./(av.^3);stdv(2) = std(ntke);
        e(2) = nanmean(ntke);
        tk = (Stat.vpro3.beam1.E(1,1:vl)+Stat.vpro3.beam3.E(1,1:vl))./2;
        av = Avgs.vpro1.x(1,1:vl)+wrms;
        ntke = tk./(av.^3);stdv(3) = std(ntke);
        e(3) = nanmean(ntke);
    end
    eps(c,:) = e;
    wvrms(c) = nanmean(wvstats.hrmsp(ind));
    stdev(c,:) = stdv;
    c = c+1;
end
eps = abs(eps);
%reorder in the order listed at the top: Q2-Q7, remember Q5 and Q6 are used
%twice (two days)
E = NaN(10,3);
E(1,:) = eps(8,:);E(2,:) = eps(6,:);E(3,:) = eps(7,:);
E(4,:) = eps(1,:);E(5,:) = eps(5,:);E(6,2) = eps(3,2);
E(7,3) = eps(3,3);E(8,2) = eps(4,2);E(9,3) = eps(4,3);
E(10,:) = eps(2,:);
WRMS = NaN(10,1);
WRMS(1) = wvrms(8);WRMS(2) = wvrms(6);WRMS(3) = wvrms(7);
WRMS(4) = wvrms(1);WRMS(5) = wvrms(5);WRMS(6) = wvrms(3);
WRMS(7) = wvrms(3);WRMS(8) = wvrms(4);WRMS(9) = wvrms(4);
WRMS(10) = wvrms(2);
STD = NaN(10,3);
STD(1,:) = stdev(8,:);STD(2,:) = stdev(6,:);STD(3,:) = stdev(7,:);
STD(4,:) = stdev(1,:);STD(5,:) = stdev(5,:);STD(6,2) = stdev(3,2);
STD(7,3) = stdev(3,3);STD(8,2) = stdev(4,2);STD(9,3) = stdev(4,3);
STD(10,:) = stdev(2,:);
%calculate veg stats
c = 1;
h = 0.005; %spacing between layers (m)
for i = 3:length(vegfiles);
    load([vegdir vegfiles{i}])
    DATA = rmfield(DATA,'info'); %remove info field for processing
    fn = fieldnames(DATA);n = length(fn);
    N = length(DATA.layer1.Rcenter);
    for j = 1:n
        [H(j),~] = size(DATA.(fn{j}).Rcenter);
        [nn,~] = size(DATA.(fn{j}).Rcenter);
        for k = 1:nn
            x = DATA.(fn{j}).Rcenter(k,1);
            y = DATA.(fn{j}).Rcenter(k,2);
            r = mean(DATA.(fn{j}).Rradii(k,:));
            lmin = DATA.(fn{j}).Layerminh;
            lmax = DATA.(fn{j}).Layermaxh;
            Vcyl(j,k) = pi*r^2*h;
            fA(j,k) = 2*r*h;
            D(j,k) = r*2;
        end
        vcyl(j,:) = sum(Vcyl(j,:));
        A(j,:) = sum(fA(j,:));
    end
    %Calculate canopy height
    H = sort(H);Hmed = find(H == round(median(H)));
    if length(Hmed) > 1
        Hmed = Hmed(ceil(end/2)); %middle median value
        hc(c) = Hmed*h; %average canopy height
    else
        hc(c) = Hmed*h; %average canopy height
    end
    %Quadrat statistics:
    phi(c) = sum(vcyl)/1; %unitless
    a(c) = sum(A)/1; %m^-1; area/1 unit volume (1m^3)
    d(c) = mean(D(1,:)); %mean diameter of layer 1
    phic(c) = (N*pi*(d(c)^2))/4; %phi modeled as a cylinder
    Ac(c) = N*d(c);
    c = c+1;
end
clear A
%recorder in the order listed at the top: Q2-Q7
Hc = NaN(10,3);
Hc(1,1:3) = hc(1);Hc(2,1:3) = hc(2);Hc(3,1:3) = hc(3);
Hc(4,1:3) = 0;Hc(5,1:3) = hc(4);Hc(6,2) = hc(5);
Hc(7,3) = hc(6);Hc(8,2) = hc(5);Hc(9,3) = hc(6);
Hc(10,1:3) = hc(7);
%vegetation descriptors
Phi = NaN(10,3);
Phi(1,1:3) = phi(1);Phi(2,1:3) = phi(2);Phi(3,1:3) = phi(3);
Phi(4,1:3) = 0;Phi(5,1:3) = phi(4);Phi(6,2) = phi(5);
Phi(7,3) = phi(6);Phi(8,2) = phi(5);Phi(9,3) = phi(6);
Phi(10,1:3) = phi(7);
Phic = NaN(10,3);
Phic(1,1:3) = phic(1);Phic(2,1:3) = phic(2);Phic(3,1:3) = phic(3);
Phic(4,1:3) = 0;Phic(5,1:3) = phic(4);Phic(6,2) = phic(5);
Phic(7,3) = phic(6);Phic(8,2) = phi(5);Phic(9,3) = phic(6);
Phic(10,1:3) = phic(7);
A = NaN(10,3);
A(1,1:3) = Ac(1);A(2,1:3) = Ac(2);A(3,1:3) = Ac(3);
A(4,1:3) = 0;A(5,1:3) = Ac(4);A(6,2) = Ac(5);
A(7,3) = Ac(6);A(8,2) = Ac(5);A(9,3) = Ac(6);
A(10,1:3) = Ac(7);
%since we're going to load the 2015 data as well, combine everything into a
%structure
dat.four.E = E;dat.four.Hc = Hc;dat.four.Phi = Phi;dat.four.Phic = Phic;
dat.four.A = A;dat.four.z = z;dat.four.wrms = WRMS;dat.four.std = STD*2;
clearvars -except dat

%%%%Load 2015 Data
veldir = 'c:\Users\bkn5\Projects\Mekong_W2015\DataAnalysis\Paper2\AveragedVelocities\';
tkedir = 'c:\Users\bkn5\Projects\Mekong_W2015\DataAnalysis\Paper2\TKE\';
vegdir = 'c:\Users\bkn5\Projects\Mekong_W2015\DataAnalysis\Paper2\VegStats\';
wvdir = 'c:\Users\bkn5\Projects\Mekong_W2015\DataAnalysis\Paper2\Spectra\';
velfiles = dir(veldir);velfiles = {velfiles.name};
tkefiles = dir(tkedir);tkefiles = {tkefiles.name};
vegfiles = dir(vegdir);vegfiles = {vegfiles.name};
wvfiles = dir(wvdir);wvfiles = {wvfiles.name};
z = [125 125 80;63 65 65;62 63 61;242 271 240;242 271 240;550 555 535;61 61 61;61 61 61;70 416 806]-43; %7 is the bin number to be used

c = 1; %index counter
eps = zeros(length(tkefiles)-2,3);
wvrms = zeros(length(tkefiles)-2,1);
stdev = zeros(length(tkefiles)-2,3);
%wvfiles are split by experiment, not by day
whichwv = [3 3 4 4 5 6 7 8 9];
for i = 3:length(tkefiles)
    load([tkedir tkefiles{i}])
    load([veldir velfiles{i}])
    k = whichwv(c);load([wvdir wvfiles{k}])
    if i == 8
        vl = 15;
    elseif i == 11
        vl = 10;
    else
        [~,vl] = size(Stat.vpro1.beam1.E);
    end
    %crop wavestats to the TKE deployment record: from closest values
    ind = zeros(1,vl);
    for j = 1:vl
        tmp = abs(wvstats.time-Stat.vpro1.time(j));
        [~,idx] = min(tmp);
        ind(j) = idx;
    end
    wvel = wvstats.rorb(ind);
    wrms = sqrt(sum(wvel.^2)/numel(wvel));
    tk = (Stat.vpro1.beam1.E(7,1:vl)+Stat.vpro1.beam3.E(7,1:vl))./2;
    av = Avgs.vpro1.x(1,1:vl)+wrms;
    ntke = tk./(av.^3);stdv(1) = std(ntke);
    e(1) = nanmean(ntke);
    tk = (Stat.vpro2.beam1.E(7,1:vl)+Stat.vpro2.beam3.E(7,1:vl))./2;
    av = Avgs.vpro2.x(1,1:vl)+wrms;
    ntke = tk./(av.^3);stdv(2) = std(ntke);
    e(2) = nanmean(ntke);
    tk = (Stat.vpro3.beam1.E(7,1:vl)+Stat.vpro3.beam3.E(7,1:vl))./2;
    av = Avgs.vpro3.x(1,1:vl)+wrms;
    ntke = tk./(av.^3);stdv(3) = std(ntke);
    e(3) = nanmean(ntke);
    eps(c,:) = e;
    wvrms(c) = nanmean(wvstats.hrmsp(ind));
    stdev(c,:) = stdv;
    c = c+1;
end
eps = abs(eps);
%reorder in the order listed at the top: Q1-Q4
E = NaN(9,3);
E(1,:) = eps(1,:);E(2,:) = eps(2,:);E(3,:) = eps(5,:);
E(4,:) = eps(6,:);E(5,:) = eps(7,:);E(6,:) = eps(8,:);
E(7,:) = eps(3,:);E(8,:) = eps(4,:);E(9,:) = eps(9,:);
WRMS = NaN(9,1);
WRMS(1) = wvrms(1);WRMS(2) = wvrms(2);WRMS(3) = wvrms(5);
WRMS(4) = wvrms(6);WRMS(5) = wvrms(7);WRMS(6) = wvrms(8);
WRMS(7) = wvrms(3);WRMS(8) = wvrms(4);WRMS(9) = wvrms(9);
STD = NaN(9,3);
STD(1,:) = stdev(1,:);STD(2,:) = stdev(2,:);STD(3,:) = stdev(5,:);
STD(4,:) = stdev(6,:);STD(5,:) = stdev(7,:);STD(6,:) = stdev(8,:);
STD(7,:) = stdev(3,:);STD(8,:) = stdev(4,:);STD(9,:) = stdev(9,:);
%calculate veg stats
c = 1;
h = 0.005; %spacing between layers (m)
for i = 3:length(vegfiles);
    load([vegdir vegfiles{i}])
    DATA = rmfield(DATA,'info'); %remove info field for processing
    fn = fieldnames(DATA);n = length(fn);
    N = length(DATA.layer1.Rcenter);
    for j = 1:n
        [H(j),~] = size(DATA.(fn{j}).Rcenter);
        [nn,~] = size(DATA.(fn{j}).Rcenter);
        for k = 1:nn
            x = DATA.(fn{j}).Rcenter(k,1);
            y = DATA.(fn{j}).Rcenter(k,2);
            r = mean(DATA.(fn{j}).Rradii(k,:));
            lmin = DATA.(fn{j}).Layerminh;
            lmax = DATA.(fn{j}).Layermaxh;
            Vcyl(j,k) = pi*r^2*h;
            fA(j,k) = 2*r*h;
            D(j,k) = r*2;
        end
        vcyl(j,:) = sum(Vcyl(j,:));
        A(j,:) = sum(fA(j,:));
    end
    %Calculate canopy height
    H = sort(H);Hmed = find(H == round(median(H)));
    if length(Hmed) > 1
        Hmed = Hmed(ceil(end/2)); %middle median value
        hc(c) = Hmed*h; %average canopy height
    else
        hc(c) = Hmed*h; %average canopy height
    end
    %Quadrat statistics:
    phi(c) = sum(vcyl)/1; %unitless
    a(c) = sum(A)/1; %m^-1; area/1 unit volume (1m^3)
    d(c) = mean(D(1,:)); %mean diameter of layer 1
    phic(c) = (N*pi*(d(c)^2))/4; %phi modeled as a cylinder
    Ac(c) = N*d(c);
    c = c+1;
end
clear A
%recorder in the order listed at the top: Q1-Q4
Hc = NaN(9,3);
Hc(1,1) = 0;Hc(1,2) = hc(1);Hc(1,3) = hc(2);
Hc(2,1) = 0;Hc(2,2) = hc(1);Hc(2,3) = hc(2);
Hc(3,:) = hc(3);Hc(4,:) = hc(3);Hc(5,:) = hc(3);
Hc(6,:) = hc(3);
Hc(7,1) = 0;Hc(7,2) = hc(4);Hc(7,3) = hc(5);
Hc(8,1) = 0;Hc(8,2) = hc(4);Hc(8,3) = hc(5);
Hc(9,:) = hc(6);
%veg statistics
Phi = NaN(9,3);
Phi(1,1) = 0;Phi(1,2) = phi(1);Phi(1,3) = phi(2);
Phi(2,1) = 0;Phi(2,2) = phi(1);Phi(2,3) = phi(2);
Phi(3,:) = phi(3);Phi(4,:) = phi(3);Phi(5,:) = phi(3);
Phi(6,:) = phi(3);
Phi(7,1) = 0;Phi(7,2) = phi(4);Phi(7,3) = phi(5);
Phi(8,1) = 0;Phi(8,2) = phi(4);Phi(8,3) = phi(5);
Phi(9,:) = phi(6);
Phic = NaN(9,3);
Phic(1,1) = 0;Phic(1,2) = phic(1);Phic(1,3) = phic(2);
Phic(2,1) = 0;Phic(2,2) = phic(1);Phic(2,3) = phic(2);
Phic(3,:) = phic(3);Phic(4,:) = phic(3);Phic(5,:) = phic(3);
Phic(6,:) = phic(3);
Phic(7,1) = 0;Phic(7,2) = phic(4);Phic(7,3) = phic(5);
Phic(8,1) = 0;Phic(8,2) = phic(4);Phic(8,3) = phic(5);
Phic(9,:) = phic(6);
A = NaN(9,3);
A(1,1) = 0;A(1,2) = Ac(1);A(1,3) = Ac(2);
A(2,1) = 0;A(2,2) = Ac(1);A(2,3) = Ac(2);
A(3,:) = Ac(3);A(4,:) = Ac(3);A(5,:) = Ac(3);
A(6,:) = Ac(3);
A(7,1) = 0;A(7,2) = Ac(4);A(7,3) = Ac(5);
A(8,1) = 0;A(8,2) = Ac(4);A(8,3) = Ac(5);
A(9,:) = Ac(6);
dat.five.E = E;dat.five.Hc = Hc;dat.five.Phi = Phi;dat.five.Phic = Phic;
dat.five.A = A;dat.five.z = z;dat.five.wrms = WRMS;dat.five.std = STD*2;
clearvars -except dat
break
%%%%END DATA CONGLOMERATOR%%%%
%Figure 1: Plot Phi against E

%Plot 2014
% zhc = (dat.four.z/1000)./dat.four.Hc; %compute z/hc, z is in mm convert to m
% zhc(4,1:3) = 0;
% %color by z/hc
% zrange = 0:0.01:2;
% color by Hs wave height
h = sqrt(2)*dat.four.wrms;
hrange = 0:0.1:0.5;
f1 = figure(1);
set(f1,'PaperOrientation','portrait',...
    'position',[400 200   600   500]);
set(gcf,'color','w','PaperPositionMode','auto')
cmap = linspecer(length(hrange),'sequential'); colormap(cmap)

for i = 1:length(dat.four.Phi)
    hold on
    for ii = 1:3
%         tmp = abs(zrange-zhc(i,ii));
        tmp = abs(hrange-h(i));
        [~,idx] = min(tmp);
        errorbar(dat.four.Phi(i,ii),log10(dat.four.E(i,ii)),...
            log10(dat.four.std(i,ii)),log10(dat.four.std(i,ii)),...
            'MarkerFaceColor',cmap(idx,:),...
            'Marker','o','MarkerSize',11,'MarkerEdgeColor','k',...
            'LineWidth',1.5,'Color','k') 
    end
end
caxis([0 0.5])
cb = colorbar;

%plot 2015
% zhc = (dat.five.z/1000)./dat.five.Hc; %compute z/hc, z is in mm convert to m
% zhc(isinf(zhc)) = 0;
h = sqrt(2)*dat.five.wrms;
hrange = 0:0.1:0.5;
cmap = linspecer(length(hrange),'sequential'); colormap(cmap)

for i = 1:length(dat.five.Phi)
    hold on
    for ii = 1:3
%         tmp = abs(zrange-zhc(i,ii));
        tmp = abs(hrange-h(i));
        [~,idx] = min(tmp);
        errorbar(dat.five.Phi(i,ii),log10(dat.five.E(i,ii)),...
            log10(dat.five.std(i,ii)),log10(dat.five.std(i,ii)),...
            'MarkerFaceColor',cmap(idx,:),...
            'Marker','o','MarkerSize',11,'MarkerEdgeColor','k',...
            'LineWidth',1.5,'Color','k') 
    end
end

%labels
xlabel('\phi','FontSize',15,'FontName','Cambria')
ylabel('log_1_0(\epsilon/(u_o_r_b+ u)^3)  (m^-^1)','FontSize',16,'FontName','Cambria')
ylabel(cb,'H_s (m)','FontSize',15,'FontName','Cambria')
title('2014-2015 Hydrodynamic Experiments','FontSize',15,'FontName','Cambria')
%global adjustments
set(gca,'LineWidth',1.5,'FontName','Cambria','FontSize',15,...
    'Xlim',[-1E-4 5E-3],'XTick',0:1E-3:5E-3,...
    'YLim',[-6 2],'YTick',-6:2:2,'box','on')
set(cb,'LineWidth',1.5,'FontName','Cambria','FontSize',15,'YTick',0:0.1:0.5)

figdir = 'c:\Users\bkn5\Projects\Mekong_W2015\Figures\Paper2\';
export_fig([figdir 'Mekong_TKEDens_cbHs'],'-jpeg','-nocrop')

%Figure 2: plot A against E
% f2 = figure(2);
% set(f2,'PaperOrientation','portrait',...
%     'position',[400 200   600   500]);
% set(gcf,'color','w','PaperPositionMode','auto')
% cmap = linspecer(length(zrange),'sequential'); colormap(cmap)
% 
% for i = 1:length(A)
%     hold on
%     for ii = 1:3
%         tmp = abs(zrange-zhc(i,ii));
%         [~,idx] = min(tmp);
%         plot(A(i,ii),log10(E(i,ii)),'MarkerFaceColor',cmap(idx,:),...
%             'Marker','o','MarkerSize',11,'MarkerEdgeColor','k',...
%             'LineWidth',1.5) 
%     end
% end
% caxis([0 2])
% cb = colorbar;
% %labels
% xlabel('Frontal Area Density (m^-^1)','FontSize',15,'FontName','Cambria')
% ylabel('log_1_0(\epsilon)  (Wkg^-^1)','FontSize',16,'FontName','Cambria')
% ylabel(cb,'z/h_c','FontSize',15,'FontName','Cambria')
% title('Autumn 2014','FontSize',15,'FontName','Cambria')
% %global adjustments
% set(gca,'LineWidth',1.5,'FontName','Cambria','FontSize',15,...
%     'Xlim',[-0.1 2],'XTick',[0:0.5:2],'box','on')
% set(cb,'LineWidth',1.5,'FontName','Cambria','FontSize',15)
% 
% % figdir = 'c:\Users\bkn5\Projects\Mekong_F2014\Figures\Paper2\';
% % export_fig([figdir 'Mekong2014_TKEDens_a_cyl'],'-jpeg','-nocrop')
% 
% %Plot Phi against Phic
% Phi(isnan(Phi)) = 0;Phic(isnan(Phic)) = 0;
% pf = polyfit(Phi(:,1),Phic(:,1),1);
% pv = polyval(pf,Phi);
% resid = Phic-pv;
% SSresid = sum(resid.^2);
% SStotal = (length(Phic-1)*var(Phic));
% rsq = 1-SSresid/SStotal;
% 
% f3 = figure(3);
% set(f3,'PaperOrientation','portrait',...
%     'position',[400 200   600   500]);
% set(gcf,'color','w','PaperPositionMode','auto')
% plot(Phi,Phic,'o','MarkerFaceColor','w',...
%     'MarkerEdgeColor','k','MarkerSize',10,...
%     'LineWidth',1.5)
% hold on
% plot(Phi,pv,'LineWidth',1.5,'Color','k')
% set(gca,'linewidth',1.5,'FontSize',12)
% xlabel('\phi')
% ylabel('\phi as cylinders')
% figdir = 'c:\Users\bkn5\Projects\Mekong_F2014\Figures\Paper2\';
% export_fig([figdir 'Mekong2014_Phi_Phic'],'-jpeg','-nocrop')
% 
% 
% 
