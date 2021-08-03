%Let's do some comparisons with other studies. For this script, I need to
%calculate Cd from our data (see Lacy et al. 2011 for method). Then, I'll
%extract Cd and Re values from other studies, and plot Cd*a by the Re #. 
clear
load('d:\Mekong_W2015\DataAnalysis\Paper2\RS_estimates\Final\RStress_10min_bdadj_f.mat')
load('D:\Mekong_W2015\DataAnalysis\Paper2\U_z_gradient_v3.mat')
load('D:\Mekong_W2015\DataAnalysis\DataReports\Vegetation\VegDat_by_height_1m.mat')

dn = fieldnames(RS);
hab = [0.062 0.063 0.061;
    0.242 0.240 0.240;
    0.550 0.550 0.550;
    0.07 0.416 0.806];
quads = {'Q2A';'Q2B';'Q2C';'Q4'};
for i = 1:4
    fn = fieldnames(RS.(dn{i}));
    for j = 1:3
        uw = RS.(dn{i}).(fn{j}).uw;
        time1 = RS.(dn{i}).(fn{j}).time;
        time2 = slope.(dn{i}).(fn{j}).time';
        U = zeros(length(time1),35);
        for kk = 1:length(time1)-1
            td = find(time2 >= time1(kk) & time2 <= time1(kk+1));
            U(kk,:) = nanmean(slope.(dn{i}).(fn{j}).u(:,td),2);
        end
        zp = hab(i,j)-0.04-linspace(0,0.03,35);
        zc = veg.five.(quads{i}).z;
        a = veg.five.(quads{i}).a;
        cid = find(zc >= zp(end) & zc <=zp(1));
        aintp = linspace(a(cid(1)),a(cid(end)),35);
        %Calculate Cd (see Lacy et al. 2011)
        %Check notebook for notation, date 04-10-17
        g = nansum(nanmean(-1*uw))-uw;
        U2 = U.^2;
        [ns,ht] = size(U2);
        gg = zeros(ns,ht);
        for kk = 1:ns
            gg(kk,:) = aintp.*U2(kk,:);
        end
        Cd = 2*g./gg;
        