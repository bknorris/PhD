%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%We're going to plot TKE dissipation rates (averaged over 10 minute
%windows) against the vegetation statistics (phi) to see if there is a
%trend in TKE with increasing vegetation density. In this case, plot PHI by
%CROSS-SHORE DISTANCE (X) and compare this to TKE dissipation rates by
%distance.

%This script is to process the vegetation data into a sensible structure.
%Distance values will be hardwired in as they are recorded from the
%deployment notes. The actual plots will be made by another script:
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%cross shore distances for Q2-Q7 (distance from the fringe). Positive X
%values are ONSHORE of the fringe, negative X values are OFFSHORE. These
%values were derived by digitizing the fringe in Google Earth.
%Data from 'QuadratDistances.xls'
X = [4.95 4.95 4.95;5.92 6.02 6.12;6.16 6.16 6.16;-37.17 -37.17 -37.17;...
    -11.48 -11.48 -11.48;10 NaN NaN;20 NaN NaN;30 NaN NaN;40 NaN NaN;...
    50 NaN NaN;60 NaN NaN;70 NaN NaN;80 NaN NaN;-10 NaN NaN;-20 NaN NaN;...
    -60.96 0 49.39;-60.96 0 49.39;11.7 11.7 11.7];
z = [99 101 98;41 43 37;83 179 279;101 204 300;84 174 272;NaN NaN NaN;
    NaN NaN NaN;NaN NaN NaN;NaN NaN NaN;NaN NaN NaN;NaN NaN NaN;NaN NaN NaN;...
    NaN NaN NaN;NaN NaN NaN;NaN NaN NaN;...
    110 42 110;110 42 110;108 208 308];
%load the 2014 Veg data
vegdir = 'd:\Projects\Mekong_F2014\DataAnalysis\Paper2\VegStats\';
vegfiles = dir([vegdir '*stats_fixed.mat']);vegfiles = {vegfiles.name};
%calculate veg stats
c = 1;
h = 0.005; %spacing between layers (m)
name = cell(length(vegfiles)-2,1);
for i = 1:length(vegfiles);
    load([vegdir vegfiles{i}])
    name{c} = regexprep(vegfiles{i},'stats_fixed.mat','');
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
    %Q5_1 and Q5_9 are 1m^2, along with the others. Q5_2-Q5_10 are 0.25m^2. 
    phi(c) = sum(vcyl)/1; %unitless
    a(c) = sum(A)/1; %m^-1; area/1 unit volume (1m^3)

    d(c) = mean(D(1,:)); %mean diameter of layer 1
    phic(c) = (N*pi*(d(c)^2))/4; %phi modeled as a cylinder
    Ac(c) = N*d(c);
    c = c+1;
end
clear A
%Organize data into a structure based on Quadrats Q2-Q7
Phi = NaN(18,3);
Phi(1,:) = phi(1);Phi(2,:) = phi(2);Phi(3,:) = phi(3);Phi(4,:) = 0;
Phi(5,:) = phi(4);Phi(6,1) = phi(6);Phi(7,1) = phi(7);Phi(8,1) = phi(8);
Phi(9,1) = phi(9);Phi(10,1) = phi(10);Phi(11,1) = phi(11);Phi(12,1) = phi(12);
Phi(13,1) = phi(13);Phi(14,1) = phi(14);Phi(15,1) = phi(5);Phi(16,1) = 0;
Phi(16,2) = phi(15);Phi(16,3) = phi(16);Phi(17,1) = 0;Phi(17,2) = phi(15);
Phi(17,3) = phi(16);
Phi(18,:) = phi(17);
Phic = NaN(18,3);
Phic(1,:) = phic(1);Phic(2,:) = phic(2);Phic(3,:) = phic(3);Phic(4,:) = 0;
Phic(5,:) = phic(4);Phic(6,1) = phic(6);Phic(7,1) = phic(7);Phic(8,1) = phic(8);
Phic(9,1) = phic(9);Phic(10,1) = phic(10);Phic(11,1) = phic(11);Phic(12,1) = phic(12);
Phic(13,1) = phic(13);Phic(14,1) = phic(14);Phic(15,1) = phic(5);Phic(16,1) = 0;
Phic(16,2) = phic(15);Phic(16,3) = phic(16);Phic(17,1) = 0;Phic(17,2) = phic(15);
Phic(17,3) = phic(16);
Phic(18,:) = phic(17);
A = NaN(18,3);
A(1,:) = a(1);A(2,:) = a(2);A(3,:) = a(3);A(4,:) = 0;
A(5,:) = a(4);A(6,1) = a(6);A(7,1) = a(7);A(8,1) = a(8);
A(9,1) = a(9);A(10,1) = a(10);A(11,1) = a(11);A(12,1) = a(12);
A(13,1) = a(13);A(14,1) = a(14);A(15,1) = a(5);A(16,1) = 0;
A(16,2) = a(15);A(16,3) = a(16);A(17,1) = 0;A(17,2) = a(15);
A(17,3) = a(16);
A(18,:) = a(17);
Hc = NaN(18,3);
Hc(1,:) = hc(1);Hc(2,:) = hc(2);Hc(3,:) = hc(3);Hc(4,:) = 0;
Hc(5,:) = hc(4);Hc(6,1) = hc(6);Hc(7,1) = hc(7);Hc(8,1) = hc(8);
Hc(9,1) = hc(9);Hc(10,1) = hc(10);Hc(11,1) = hc(11);Hc(12,1) = hc(12);
Hc(13,1) = hc(13);Hc(14,1) = hc(14);Hc(15,1) = hc(5);Hc(16,1) = 0;
Hc(16,2) = hc(15);Hc(16,3) = hc(16);Hc(17,1) = 0;Hc(17,2) = hc(15);
Hc(17,3) = hc(16);
Hc(18,:) = hc(17);
mud = {'Mudflat'};
quad = cell(18,3);
quad(1,:) = name(1);quad(2,:) = name(2);quad(3,:) = name(3);quad(4,:) = mud;
quad(5,:) = name(4);quad(6,1) = name(6);quad(7,1) = name(7);quad(8,1) = name(8);
quad(9,1) = name(9);quad(10,1) = name(10);quad(11,1) = name(11);quad(12,1) = name(12);
quad(13,1) = name(13);quad(14,1) = name(14);quad(15,1) = name(5);quad(16,1) = mud;
quad(16,2) = name(15);quad(16,3) = name(16);quad(17,1) = mud;quad(17,2) = name(15);
quad(17,3) = name(16);
quad(18,:) = name(17);
vegdat.four.quad = quad;vegdat.four.hc = Hc;
vegdat.four.Phi = Phi;vegdat.four.Phic = Phic;
vegdat.four.A = A;vegdat.four.X = X;vegdat.four.Z = z;
clearvars -except vegdat veldat h

%%%%Load 2015 Data
X = [-50.54 0.54 93.76;-50.54 0.54 93.76;0 0 0;0 0 0;0 0 0;0 0 0;...
    -50.35 1.37 72.64;-50.35 1.37 72.64;24.52 24.52 24.52];
z = [125 125 80;63 65 65;62 63 61;242 271 240;242 271 240;...
    550 555 535;61 61 61;61 61 61;70 416 806]-43; %7 is the bin number to be used
vegdir = 'd:\Projects\Mekong_W2015\DataAnalysis\Paper2\VegStats\';
vegfiles = dir([vegdir '*stats_fixed.mat']);vegfiles = {vegfiles.name};
%calculate veg stats
c = 1;
for i = 1:length(vegfiles);
    load([vegdir vegfiles{i}])
    DATA = rmfield(DATA,'info'); %remove info field for processing
    name{c} = regexprep(vegfiles{i},'stats_fixed.mat','');
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
%Organize data into a structure based on Quadrats Q2-Q7
Hc = NaN(9,3);
Hc(1,1) = 0;Hc(1,2) = hc(1);Hc(1,3) = hc(2);
Hc(2,1) = 0;Hc(2,2) = hc(1);Hc(2,3) = hc(2);
Hc(3,:) = hc(3);Hc(4,:) = hc(4);Hc(5,:) = hc(5);
Hc(6,:) = hc(5);
Hc(7,1) = 0;Hc(7,2) = hc(6);Hc(7,3) = hc(7);
Hc(8,1) = 0;Hc(8,2) = hc(6);Hc(8,3) = hc(7);
Hc(9,:) = hc(8);
%veg statistics
Phi = NaN(9,3);
Phi(1,1) = 0;Phi(1,2) = phi(1);Phi(1,3) = phi(2);
Phi(2,1) = 0;Phi(2,2) = phi(1);Phi(2,3) = phi(2);
Phi(3,:) = phi(3);Phi(4,:) = phi(4);Phi(5,:) = phi(5);
Phi(6,:) = phi(5);
Phi(7,1) = 0;Phi(7,2) = phi(6);Phi(7,3) = phi(7);
Phi(8,1) = 0;Phi(8,2) = phi(6);Phi(8,3) = phi(7);
Phi(9,:) = phi(8);
Phic = NaN(9,3);
Phic(1,1) = 0;Phic(1,2) = phic(1);Phic(1,3) = phic(2);
Phic(2,1) = 0;Phic(2,2) = phic(1);Phic(2,3) = phic(2);
Phic(3,:) = phic(3);Phic(4,:) = phic(4);Phic(5,:) = phic(5);
Phic(6,:) = phic(5);
Phic(7,1) = 0;Phic(7,2) = phic(6);Phic(7,3) = phic(7);
Phic(8,1) = 0;Phic(8,2) = phic(6);Phic(8,3) = phic(7);
Phic(9,:) = phic(8);
A = NaN(9,3);
A(1,1) = 0;A(1,2) = a(1);A(1,3) = a(2);
A(2,1) = 0;A(2,2) = a(1);A(2,3) = a(2);
A(3,:) = a(3);A(4,:) = a(4);A(5,:) = a(5);
A(6,:) = a(5);
A(7,1) = 0;A(7,2) = a(6);A(7,3) = a(7);
A(8,1) = 0;A(8,2) = a(6);A(8,3) = a(7);
A(9,:) = a(8);
mud = {'Mudflat'};
quad = cell(9,3);
quad(1,1) = mud;quad(1,2) = name(1);quad(1,3) = name(2);
quad(2,1) = mud;quad(2,2) = name(1);quad(2,3) = name(2);
quad(3,:) = name(3);quad(4,:) = name(4);quad(5,:) = name(5);
quad(6,:) = name(5);
quad(7,1) = mud;quad(7,2) = name(6);quad(7,3) = name(7);
quad(8,1) = mud;quad(8,2) = name(6);quad(8,3) = name(7);
quad(9,:) = name(8);

vegdat.five.quad = quad;
vegdat.five.Hc = Hc;vegdat.five.Phi = Phi;
vegdat.five.Phic = Phic;vegdat.five.A = A;
vegdat.five.X = X;vegdat.five.Z = z;
clearvars -except vegdat veldat 
save(['d:\Projects\Mekong_W2015\DataAnalysis\Paper2\Environment\' 'VegetationData'],'vegdat','-v7.3')
